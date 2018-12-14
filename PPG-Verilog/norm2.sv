`timescale 1ns / 1ns

module norm2#(parameter I = 20,
			  parameter J = 240,
			  parameter Q = 15,
			  parameter N = 32)(
			  input clk,
			  input rst_n,
			  input[N-1:0] A[0:I-1][0:J-1],
			  input start,
			  output[N-1:0] A_norm2[0:J-1],
			  output done);

wire[J-1:0] complete; 			  
genvar j, c, d;
wire[N-1:0] col_A[0:J-1][0:I-1];

generate
// Generate col_A, columns of A are rows of A_in
	for(c = 0; c < J; c = c + 1) begin: GEN_COL_A_0
		for(d = 0; d < I; d = d + 1) begin: GEN_COL_A_1
			assign col_A[c][d] = A[d][c];
		end
	end

// Instantiate J Square_Sum Units, one for each element of A_norm2
	for(j = 0; j < J; j = j+1) begin: A_COLUMNS
		Square_Sum #(.I(I),
					 .Q(Q),
					 .N(N))iSS(.clk(clk),
							  .rst_n(rst_n),
							  .start(start),
							  .in(col_A[j]),
							  .out(A_norm2[j]),
							  .done(complete[j]));
	end
endgenerate

assign done = &complete;										// done asserted when all units are complete

endmodule