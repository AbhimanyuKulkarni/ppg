`timescale 1ns / 1ns

module dot#(parameter I = 20,							
			parameter J = 240,
			parameter Q = 15,
			parameter N = 32)(
			input[N-1:0] in_A[0:I-1][0:J-1], 				// in_A --> IxJ 
			input[N-1:0] in_B[0:J-1],						// in_B --> Jx1 
			input start,
			input clk,
			input rst_n,
			output done,
			output[N-1:0] out_C[0:I-1]);					// out_C --> Ix1

reg[0:I-1] complete;
	
genvar i;

generate
// Instantiate I MAC units, one for each element of out_C
	for(i = 0; i < I; i = i + 1) begin: A_ROWS				// i traverses over rows of in_A
			MAC #(.Q(Q),
				  .N(N),
				  .J(J))iMAC(.in_A(in_A[i]), 				// i'th row of in_A
						   .in_B(in_B),					
						   .clk(clk),
						   .rst_n(rst_n),
						   .start(start),
						   .out_C(out_C[i]),
						   .done(complete[i]));	
	end
endgenerate

assign done = &complete;									// done is asserted when complete is all 1s, i.e. all MAC units are done.
	 
endmodule