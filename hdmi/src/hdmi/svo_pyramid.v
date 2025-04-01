`timescale 1ns / 1ps
`include "hdmi/svo_defines.vh"

module svo_pyramid #(
    // Pass the same parameters as in svo_tcard or from your SVO system
    parameter SVO_MODE             = "640x480V",
    parameter SVO_FRAMERATE        = 60,
    parameter SVO_BITS_PER_PIXEL   = 24,
    parameter SVO_BITS_PER_RED     = 8,
    parameter SVO_BITS_PER_GREEN   = 8,
    parameter SVO_BITS_PER_BLUE    = 8,
    parameter SVO_BITS_PER_ALPHA   = 0
)(
    input  wire clk,
    input  wire resetn,

    // AXI-Stream output interface
    output reg                    out_axis_tvalid,
    input  wire                   out_axis_tready,
    output reg [SVO_BITS_PER_PIXEL-1:0] out_axis_tdata,
    output reg [0:0]             out_axis_tuser
);

    //---------------------------------------------------------
    // 1) Basic timing/counters to produce 640×480 pixels
    //---------------------------------------------------------
    // For simplicity, assume we produce a pixel every clock
    // for the active region only, ignoring blanking intervals.
    //
    // In a real design, you might match the horizontal/vertical
    // timing from SVO or replicate the behavior of svo_tcard.
    //---------------------------------------------------------

    localparam H_RES = 640;
    localparam V_RES = 480;

    reg [9:0] x;  // Enough bits for 640
    reg [8:0] y;  // Enough bits for 480

    // 'angle' increments each frame to spin the pyramid
    reg [15:0] angle;

    // When we finish the last pixel of a frame, increment 'angle'
    wire last_pixel = (x == H_RES-1) && (y == V_RES-1);

    always @(posedge clk) begin
        if (!resetn) begin
            x <= 0;
            y <= 0;
            angle <= 0;
        end else begin
            if (out_axis_tready) begin
                // Move to next pixel
                if (x == H_RES-1) begin
                    x <= 0;
                    if (y == V_RES-1) begin
                        y <= 0;
                        angle <= angle + 32;  // spin once per frame
                    end else begin
                        y <= y + 1;
                    end
                end else begin
                    x <= x + 1;
                end
            end
        end
    end

    //---------------------------------------------------------
    // 2) 3D Tetrahedron geometry
    //---------------------------------------------------------
    // We'll define four vertices in "object space". For example,
    // a simple tetrahedron (pyramid) centered near the origin:
    //
    // v0 = (  0,  50,   0 )
    // v1 = ( 50, -50,  50 )
    // v2 = (-50, -50,  50 )
    // v3 = (  0, -50, -50 )
    //
    // We'll rotate them around the Y-axis by 'angle' each frame,
    // then apply a simple perspective projection onto 2D.
    //---------------------------------------------------------

    // For a small example, let's do integer or simple fixed point.
    // We'll define a little sine/cosine table or we can do a
    // minimal approximate approach.

    // A quick approximate LUT for sin/cos by angle[7:0], ignoring
    // fractional bits for brevity. (In a real design, you'd want
    // more precision, or a cordic, etc.)
    reg signed [15:0] sinLUT [0:255];
    reg signed [15:0] cosLUT [0:255];

    initial begin : init_luts
        integer i;
        for (i = 0; i < 256; i = i+1) begin
            // Multiply by 128 or so to keep some fractional bits
            sinLUT[i] = $rtoi(128.0 * $sin(2.0 * 3.14159 * i / 256.0));
            cosLUT[i] = $rtoi(128.0 * $cos(2.0 * 3.14159 * i / 256.0));
        end
    end

    // Our base vertices (object space)
    localparam signed [15:0] V0X =  0,  V0Y =  50, V0Z =   0;
    localparam signed [15:0] V1X =  50, V1Y = -50, V1Z =  50;
    localparam signed [15:0] V2X = -50, V2Y = -50, V2Z =  50;
    localparam signed [15:0] V3X =  0,  V3Y = -50, V3Z = -50;

    //---------------------------------------------------------
    // 3) Projection pipeline (wireframe)
    //---------------------------------------------------------
    // On every pixel, we do:
    //   1) Rotate each vertex by 'angle'.
    //   2) Perspective project to 2D.
    //   3) Check if the pixel is near any of the edges.
    //---------------------------------------------------------

    // We'll do everything combinational for illustration.
    // For a real design, you'd likely pipeline these steps.

    // Center of screen = (H_RES/2, V_RES/2)
    localparam signed [10:0] X_CENTER = H_RES/2;
    localparam signed [9:0]  Y_CENTER = V_RES/2;

    // Simple perspective constants
    localparam signed [15:0] FOCAL_LEN = 128;  // bigger -> less perspective

    // Return a small integer distance from (px,py) to the line
    // from (x1,y1) to (x2,y2).
    function [15:0] line_dist;
        input signed [10:0] px, py;
        input signed [10:0] x1, y1;
        input signed [10:0] x2, y2;
        // We do a minimal integer distance = area*2 / length.
        reg signed [31:0] dx, dy;
        reg signed [31:0] num, den;
    begin
        dx = x2 - x1;
        dy = y2 - y1;
        num = (px - x1)*dy - (py - y1)*dx;  // cross product
        if (num < 0) num = -num;           // absolute value
        den = (dx*dx + dy*dy);
        // avoid divide-by-zero
        if (den == 0) begin
            line_dist = 16'hFFFF;
        end else begin
            // approximate integer version of area / sqrt(den)
            // do a rough scale to keep it small
            // area = num, length = sqrt(den) ~ avoid real sqrt for example
            // We'll just do area^2 / den as a rough measure
            line_dist = num[15:0];  // very rough or do better
        end
    end
    endfunction

    // For quick reading of sin/cos
    wire [7:0]  angle_idx = angle[15:8];  // top bits
    wire signed [15:0] s = sinLUT[ angle_idx ];
    wire signed [15:0] c = cosLUT[ angle_idx ];

    // Let’s define a small task that:
    //   - rotates (vx,vy,vz) about Y
    //   - then does perspective project
    //   - returns screen coords (sx, sy)
    task project_vertex(
        input  signed [15:0] vx,
        input  signed [15:0] vy,
        input  signed [15:0] vz,
        output signed [10:0] sx,
        output signed [9:0]  sy
    );
        reg signed [31:0] rx, ry, rz;  // rotated coords
        reg signed [31:0] px, py;      // projected coords
        reg signed [31:0] denom;
        reg signed [31:0] scale_Q8; // stored as unsigned bits
        begin
            // rotate about Y
            //  rx = c*vx + 0*vy + -s*vz
            //  ry = 1*vy   (no rotation in Y about Y axis)
            //  rz = s*vx + c*vz
            rx = (c * vx - s * vz) >>> 7;  // shift down by 7 because s,c scaled by 128
            ry = vy; // no change in Y for Y rotation
            rz = (s * vx + c * vz) >>> 7;

            // perspective scale = FOCAL_LEN / (FOCAL_LEN - rz)
            // We'll do an integer approximation:
            // screen_x = X_CENTER + rx * scale
            // screen_y = Y_CENTER - ry * scale   (minus because +Y down in screen coords)
            // scale = FOCAL_LEN/(FOCAL_LEN - rz)
            // We'll do fixed-point: scale_Q8 = (FOCAL_LEN << 8)/(FOCAL_LEN - rz)
          

            denom = (FOCAL_LEN - rz);
            if (denom == 0) denom = 1; // avoid /0
            scale_Q8 = ((FOCAL_LEN << 8) / denom);

            px = X_CENTER + ((rx * scale_Q8) >>> 8);
            py = Y_CENTER - ((ry * scale_Q8) >>> 8);

            sx = px[10:0];
            sy = py[9:0];
        end
    endtask

    // The edges of the tetrahedron
    // (v0-v1, v0-v2, v0-v3, v1-v2, v2-v3, v3-v1)
    reg signed [10:0] sx0, sx1, sx2, sx3; 
    reg signed [9:0]  sy0, sy1, sy2, sy3;

    // We'll update the projected vertices any time x=0,y=0 (start of frame),
    // so they're stable during the entire frame. That means the pyramid
    // doesn’t keep rotating per-pixel, only once per frame.
    always @* begin
        // Default to old values, but if we’re at the start of frame,
        // recalc them. For simplicity, we'll do it combinationally,
        // but you could do a small state machine, etc.
        sx0 = 0; sy0 = 0;
        sx1 = 0; sy1 = 0;
        sx2 = 0; sy2 = 0;
        sx3 = 0; sy3 = 0;
    end

    // A small register to hold them for the entire frame
    reg signed [10:0] sx0_r, sx1_r, sx2_r, sx3_r; 
    reg signed [9:0]  sy0_r, sy1_r, sy2_r, sy3_r;

    always @(posedge clk) begin
        if (!resetn) begin
            sx0_r <= 0; sy0_r <= 0;
            sx1_r <= 0; sy1_r <= 0;
            sx2_r <= 0; sy2_r <= 0;
            sx3_r <= 0; sy3_r <= 0;
        end else if (x == 0 && y == 0 && out_axis_tready) begin
            // Recompute once per frame
            project_vertex(V0X, V0Y, V0Z, sx0_r, sy0_r);
            project_vertex(V1X, V1Y, V1Z, sx1_r, sy1_r);
            project_vertex(V2X, V2Y, V2Z, sx2_r, sy2_r);
            project_vertex(V3X, V3Y, V3Z, sx3_r, sy3_r);
        end
    end

    //---------------------------------------------------------
    // 4) Decide pixel color
    //---------------------------------------------------------
    // We'll do a quick “distance to each edge.” If < threshold,
    // draw white; else black.
    //---------------------------------------------------------

    // Edges: (0-1), (0-2), (0-3), (1-2), (2-3), (3-1)
    // We'll check them all and keep the minimum distance.
    wire [15:0] dist01, dist02, dist03, dist12, dist23, dist31;
    assign dist01 = line_dist(x, y, sx0_r, sy0_r, sx1_r, sy1_r);
    assign dist02 = line_dist(x, y, sx0_r, sy0_r, sx2_r, sy2_r);
    assign dist03 = line_dist(x, y, sx0_r, sy0_r, sx3_r, sy3_r);
    assign dist12 = line_dist(x, y, sx1_r, sy1_r, sx2_r, sy2_r);
    assign dist23 = line_dist(x, y, sx2_r, sy2_r, sx3_r, sy3_r);
    assign dist31 = line_dist(x, y, sx3_r, sy3_r, sx1_r, sy1_r);

    wire [15:0] min_d1 = (dist01 < dist02) ? dist01 : dist02;
    wire [15:0] min_d2 = (dist03 < dist12) ? dist03 : dist12;
    wire [15:0] min_d3 = (dist23 < dist31) ? dist23 : dist31;
    wire [15:0] min_d4 = (min_d1 < min_d2) ? min_d1 : min_d2;
    wire [15:0] min_dist = (min_d3 < min_d4) ? min_d3 : min_d4;

    // threshold for “on edge”
    localparam [15:0] EDGE_THRESHOLD = 200; // tune as needed

    // If within threshold, color = white. Otherwise black.
    // We'll form a 24-bit color: {8'dR, 8'dG, 8'dB}.
    wire [7:0] color_r = (min_dist < EDGE_THRESHOLD) ? 8'hFF : 8'h00;
    wire [7:0] color_g = (min_dist < EDGE_THRESHOLD) ? 8'hFF : 8'h00;
    wire [7:0] color_b = (min_dist < EDGE_THRESHOLD) ? 8'hFF : 8'h00;

    // Combine into out_axis_tdata
    wire [SVO_BITS_PER_PIXEL-1:0] pixel_color = { color_r, color_g, color_b };

    //---------------------------------------------------------
    // 5) Drive the AXI-Stream outputs
    //---------------------------------------------------------

    always @(posedge clk) begin
        if (!resetn) begin
            out_axis_tvalid <= 1'b0;
            out_axis_tdata  <= {SVO_BITS_PER_PIXEL{1'b0}};
            out_axis_tuser  <= 1'b0;
        end else begin
            // We produce valid data whenever x < H_RES and y < V_RES.
            // Once we've reached the last pixel, we remain valid,
            // or we could set it to 0 in blanking. 
            // For simplicity, let's just keep it valid for the entire
            // active region. Real designs might also handle blanking times.
            out_axis_tvalid <= (x < H_RES && y < V_RES);

            // Mark start-of-frame (TUSER=1) on the very first pixel.
            // Otherwise 0.
            if (x == 0 && y == 0) 
                out_axis_tuser <= 1'b1;
            else
                out_axis_tuser <= 1'b0;

            if (out_axis_tready && out_axis_tvalid) begin
                out_axis_tdata <= pixel_color;
            end
        end
    end

endmodule