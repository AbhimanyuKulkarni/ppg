`timescale 1ns / 1ns

module calc_nxt_xhat#(parameter Q = 15,
					  parameter N = 32)(
					  input clk,
					  input rst_n,
					  input start,
					  input[N-1:0] rho,
					  input[N-1:0] lambda,
					  input[N-1:0] A_norm2,
					  output reg[N-1:0] nxt_xhat,
					  output ready);
					  
reg activate;											// To activate the qdiv module

wire[N-1:0] dividend;
wire[N-1:0] diff;										// |rho| - lambda
wire[N-1:0] neg_lambda;									// Negative of lambda
wire overflow;
reg sign;												// Store sign of rho separately
reg[N-1:0] rho_int;										// Sample the value of rho at start
wire[N-1:0] new_xhat;

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		rho_int <= 1'b0;
	else if(start)
		rho_int <= rho;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		nxt_xhat <= 'h0;
	else if(start)
		nxt_xhat <= 'h0;
	else
		nxt_xhat <= new_xhat;
end

always_ff@(posedge clk) begin
	sign <= rho[N-1];
	rho_int[N-1] <= 1'b0;									// Make rho positive (take absolute value)
end

assign neg_lambda[N-2:0] = lambda[N-2:0];				// Copy the magnitude
assign neg_lambda[N-1] = ~lambda[N-1];					// Invert the sign

qadd #(.Q(Q),
	   .N(N))iSUB(.a(rho_int),
				  .b(neg_lambda), 
				  .c(diff));
				  
assign dividend[N-2:0] = (diff[N-1]) ? 0 : diff[N-2:0];	// Magnitude of dividend is diff if diff is positive, otherwise 0
assign dividend[N-1] = sign ? ~diff[N-1] : diff[N-1];	// Invert the sign of diff if rho is negative, otherwise keep the same
				  
always_ff@(posedge clk or negedge rst_n) begin			// Division starts one clock cycle after start input is asserted to
	if(!rst_n)											// give enough time for all inputs to be prepared.
		activate <= 1'b0;
	else
		activate <= start;
end

qdiv #(.Q(Q),
	   .N(N))iDIV(.i_clk(clk),
				  .rst_n(rst_n),
				  .i_start(activate),
				  .i_dividend(dividend),
				  .i_divisor(A_norm2),
				  .o_quotient_out(new_xhat),
				  .o_complete(ready),
				  .o_overflow(overflow),
				  .div_0());

endmodule