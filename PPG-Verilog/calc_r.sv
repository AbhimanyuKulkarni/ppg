`timescale 1ns / 1ns

module calc_r#(parameter I = 20,
			   parameter Q = 15,
			   parameter N = 32)(
			   input clk,
			   input rst_n,
			   input start,
			   input op,										// op=0 : nxt_r = r - A*x_hat
			   input[N-1:0] A[0:I-1],							// op=1 : nxt_r = r + A*x_hat
			   input[N-1:0] r[0:I-1],
			   input[N-1:0] x_hat,
			   output reg[N-1:0] nxt_r[0:I-1],
			   output reg ready);

wire[I-1:0] complete;
wire[I-1:0] overflows;
wire[N-1:0] dr[0:I-1];
wire[N-1:0] dr_add[0:I-1]; 
wire[N-1:0] new_r[0:I-1];
reg[N-1:0] r_in[0:I-1];

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r_in <= '{default:'h0};
	else if(start)
		r_in <= r;
end
			   
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		nxt_r <= '{default:'h0};
	else
		nxt_r <= new_r;
end
			   
genvar i;

generate														// Instantiate I qmults and qadd units, one for each element of nxt_r
	for(i = 0; i < I; i = i + 1) begin: QMULTS
		qmults #(.Q(Q),
				 .N(N))iQMULTS(.i_clk(clk),
							  .rst_n(rst_n),
							  .i_multiplicand(A[i]),
							  .i_multiplier(x_hat),
							  .i_start(start),
							  .o_result_out(dr[i]),
							  .o_complete(complete[i]),
							  .o_overflow(overflows[i]));
		
		assign dr_add[i][N-1] = op ? dr[i][N-1] : ~dr[i][N-1];		// Negate dr if op = 0
		assign dr_add[i][N-2:0] = dr[i][N-2:0];
		
		qadd #(.Q(Q),
			   .N(N))iQADD(.a(r_in[i]), 
						  .b(dr_add[i]), 
						  .c(new_r[i]));
	end						  
endgenerate

always_ff@(posedge clk or negedge rst_n) begin					// ready is flopped to give time for qadd to complete their operation
	if(!rst_n)
		ready <= 1'b0;
	else
		ready <= &complete;
end

endmodule
					