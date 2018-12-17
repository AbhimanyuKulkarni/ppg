`timescale 1ns / 1ns

module vecsub#(parameter I = 20,
			   parameter Q = 15,
			   parameter N = 32)(
			   input[N-1:0] in1[0:I-1],
			   input[N-1:0] in2[0:I-1],
			   output[N-1:0] out[0:I-1]);

wire[N-1:0] in2_neg[0:I-1];									// The negative of in2, which we will add to in1
			   			   
genvar i;

generate													// Instantiate I qadd modules, one for each element of out 
	for(i = 0; i < I; i = i+1) begin: SUB
		assign in2_neg[i][N-2:0] = in2[i][N-2:0];			// Copy the value of in2 to in2_neg
		assign in2_neg[i][N-1] = ~in2[i][N-1];				// Negate the sign bit of in2
		qadd #(.Q(Q),
			   .N(N))iqadd(.a(in1[i]),
						  .b(in2_neg[i]),
						  .c(out[i]));
	end
endgenerate

endmodule