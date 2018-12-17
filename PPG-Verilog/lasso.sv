`timescale 1ns / 1ns

module lasso#(parameter I = 20,
			  parameter J = 240,
			  parameter Q = 7,
			  parameter N = 16)(
			  input clk,
			  input rst_n,
			  input start,
			  input[N-1:0] y[0:I-1],
			  input[N-1:0] A[0:I-1][0:J-1],
			  input[N-1:0] lambda,
			  input[N-1:0] tol,
			  output reg[N-1:0] xhat[0:J-1],
			  output reg done);

reg init, renew, iterate, div, dot;
wire norm2_done, iterate_done, div_done, dot_done1;
reg dot_done;
reg set_done, clr_done;
wire[N-1:0] A_norm2[0:J-1];
reg[N-1:0] max_xj, max_dxj;
wire[N-1:0] max_xj_new, max_dxj_new;
reg[N-1:0] r[0:I-1];
wire[N-1:0] r_init[0:I-1];
wire[N-1:0] r_new[0:I-1];
wire[N-1:0] xhat_new[0:J-1];
wire[N-1:0] quotient;
wire div_0;
wire overflow;
//wire nxt_stop;
reg stop;
wire[N-1:0] yhat[0:I-1];

typedef enum logic[2:0]{IDLE, DOT, INIT, ITERATE, RENEW, DIV, CHECK} state_t;
state_t state, nxt_state; 
			  
// Initialise norm2 for calculating A_norm2
norm2 #(.I(I),
		.J(J),
		.Q(Q),
		.N(N)) iAnorm2(.clk(clk),
					   .rst_n(rst_n),
					   .start(init),
					   .A(A),
					   .A_norm2(A_norm2),
					   .done(norm2_done));
					   
// Holding registers for max_xj, max_dxj, r and xhat
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_xj <= 0;
	else if(init)
		max_xj <= 0;
	else if(renew)
		max_xj <= max_xj_new;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_dxj <= 0;
	else if(init)
		max_dxj <= 0;
	else if(renew)
		max_dxj <= max_dxj_new;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r <= '{default:0};
	else if(init)
		r <= r_init;														// Initial xhat = 0, yhat = 0, r = y-yhat = y.
	else if(renew)
		r <= r_new;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		xhat <= '{default:0};
	else if(start)
		xhat <= '{default:4};
	else if(renew)
		xhat <= xhat_new;
end

// Calculate yhat
dot #(.I(I),
	  .J(J),
	  .Q(Q),
	  .N(N)) idot(.clk(clk),
				  .rst_n(rst_n),
				  .start(dot),
				  .in_A(A),
				  .in_B(xhat),
				  .out_C(yhat),
				  .done(dot_done1));
				  
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		dot_done <= 1'b0;
	else
		dot_done <= dot_done1;
end

// Calculate initial r
vecsub #(.I(I),
		 .Q(Q),
		 .N(N)) ivecsub(.in1(y),
						.in2(yhat),
						.out(r_init));				  
				  
// Initialise Iterate for doing the iterations
Iterate #(.I(I),
		  .J(J),
		  .Q(Q),
		  .N(N)) iIterate(.clk(clk),
						  .rst_n(rst_n),
						  .start(iterate),
						  .xhat_in(xhat),
						  .A(A),
						  .A_norm2_in(A_norm2),
						  .max_xj_in(max_xj),
						  .r_in(r),
						  .lambda(lambda),
						  .r_out(r_new),
						  .max_xj_out(max_xj_new),
						  .max_dxj_out(max_dxj_new),
						  .xhat_out(xhat_new),
						  .done(iterate_done));

// Initialise qdiv for calculating max_dxj/max_xj
qdiv #(.Q(Q),
	   .N(N)) iqdiv(.i_dividend(max_dxj),
					.i_divisor(max_xj),
					.i_start(div),
					.i_clk(clk),
					.rst_n(rst_n),
					.o_quotient_out(quotient),
					.o_complete(div_done),
					.o_overflow(overflow),
					.div_0(div_0));

/*always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		stop <= 1'b0;
	else if(div_done)
		stop <= nxt_stop;
end*/
					
max #(.Q(Q), 
	  .N(N)) imax(.a(quotient),
				  .b(tol),
				  .max(),
				  .which(stop));						// If tol is greater than quotient, stop = 1

// Infer the state flop
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		state <= IDLE;
	else 
		state <= nxt_state;
end

// Next state and output logic
always_comb begin
	nxt_state = IDLE;
	dot = 1'b0;
	init = 1'b0;
	iterate = 1'b0;
	renew = 1'b0;
	div = 1'b0;
	set_done = 1'b0;
	clr_done = 1'b0;
	
	case(state)
		IDLE:	begin
					if(start) begin
						dot = 1'b1;
						clr_done = 1'b1;
						nxt_state = DOT;
					end
					else
						nxt_state = IDLE;
				end
				
		DOT:	begin
					if(dot_done) begin
						init = 1'b1;
						nxt_state = INIT;
					end
					else
						nxt_state = DOT;
				end
				
		INIT: 	begin
					if(norm2_done) begin
						iterate = 1'b1;
						nxt_state = ITERATE;
					end
					else
						nxt_state = INIT;
				end
				
		ITERATE:begin
					if(iterate_done) begin
						renew = 1'b1;
						nxt_state = RENEW;
					end
					else
						nxt_state = ITERATE;
				end
				
		RENEW:	begin
					div = 1'b1;
					nxt_state = DIV;
				end
		
		DIV: 	begin
					if(div_done)
						nxt_state = CHECK;
					else
						nxt_state = DIV;
				end
			
		CHECK:	begin
					if(stop) begin
						if(div_0) begin
							iterate = 1'b1;
							nxt_state = ITERATE;
						end
						else begin
							set_done = 1'b1;
							nxt_state = IDLE;
						end
					end
					else begin
						iterate = 1'b1;
						nxt_state = ITERATE;
					end
				end
				
		default:	nxt_state = IDLE;
						
	endcase
end

// S-R Flip Flop for done
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		done <= 1'b0;
	else if(clr_done)
		done <= 1'b0;
	else if(set_done)
		done <= 1'b1;
end

endmodule	