`timescale 1ns / 1ns

module Update_Element#(parameter I = 20,
					   parameter Q = 15,
					   parameter N = 32)(
					   input clk,
					   input rst_n,
					   input start,
					   input[N-1:0] xhat_j,
					   input[N-1:0] A_norm2_j,
					   input[N-1:0] lambda,
					   input[N-1:0] A_j[0:I-1],
					   input[N-1:0] max_xj_in,
					   input[N-1:0] max_dxj_in,
					   input[N-1:0] r_in[0:I-1],
					   output reg[N-1:0] r_out[0:I-1],
					   output reg[N-1:0] max_dxj_out,
					   output reg[N-1:0] max_xj_out,
					   output reg[N-1:0] nxt_xhat_j,
					   output reg done);
					  
typedef enum logic[2:0] {IDLE, CHECK, CALC_R_INT, CALC_RHO, CALC_X, CALC_R_OUT, WAIT} state_t;
state_t state, nxt_state;

reg start_r_int;
wire r_int_done;
//wire[N-1:0] nxt_r_int[0:I-1];									// To store the intermediate value of r
reg[N-1:0] r_int[0:I-1];
//wire[N-1:0] nxt_r_out[0:I-1];
reg[N-1:0] rho_j;// nxt_rho_j;
wire rho_done;
reg start_rho;
wire[N-1:0] new_xhat_j;
//reg[N-1:0] r_input[0:I-1];
//reg[N-1:0] A_input[0:I-1];
//reg[N-1:0] A_norm2_input;

reg rst_init;
reg calc_max, calc_max1, update_max;

reg start_xhat;
wire xhat_done;

reg start_r_nxt;
wire r_nxt_done;

wire[N-1:0] neg_xhat_j;
wire[N-1:0] dxj;
//wire[N-1:0] nxt_max_xj_out, nxt_max_dxj_out;
//reg[N-1:0] a_xj, b_xj, a_dxj, b_dxj;

reg set_done, clr_done;

////// Holding registers //////
/*always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r_input <= '{default:'h0};
	else if(start)
		r_input <= r_in;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		A_input <= '{default:'h0};
	else if(start)
		A_input <= A_j;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		A_norm2_input <= 'h0;
	else if(start)
		A_norm2_input <= A_norm2_j;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		rho_j <= 'h0;
	else if(rho_done)
		rho_j <= nxt_rho_j;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r_int <= '{default:'h0};
	else if(r_int_done)
		r_int <= nxt_r_int;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r_out <= '{default:'h0};
	else if(r_nxt_done)
		r_out <= nxt_r_out;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		nxt_xhat_j <= 'h0;
	else if(start)
		nxt_xhat_j <= 'h0;
	else if(xhat_done)
		nxt_xhat_j <= new_xhat_j;
end*/


////// Instantiate calc_r module to get intermediate value of r ///////
calc_r #(.I(I), .Q(Q), .N(N)) icalc_int_r(.clk(clk),
										  .rst_n(rst_n),
										  .start(start_r_int),
										  .op(1'b1),
										  .A(A_j),
										  .r(r_in),
										  .x_hat(xhat_j),
										  .nxt_r(r_int),
										  .ready(r_int_done));				 
									  
////// Instantiate MAC for calculation of rho_j. Parameter J is given the value I //////
MAC #(.J(I), .Q(Q), .N(N)) iMAC_rho(.clk(clk),
								    .rst_n(rst_init),
									.start(start_rho),
									.in_A(A_j),
									.in_B(r_int),
									.out_C(rho_j),
									.done(rho_done));
									  
////// Instantiate calc_nxt_xhat to calculate new value of xhat ///////
calc_nxt_xhat #(.Q(Q), .N(N)) icalc_xhat(.clk(clk),
										 .rst_n(rst_init),
										 .start(start_xhat),
										 .rho(rho_j),
										 .lambda(lambda),
										 .A_norm2(A_norm2_j),
										 .nxt_xhat(nxt_xhat_j),
										 .ready(xhat_done));
										 
/////// Instantiate calc_r again to compute the new value of r ///////
calc_r #(.I(I), .Q(Q), .N(N)) icalc_nxt_r(.clk(clk),
										  .rst_n(rst_n),
										  .start(start_r_nxt),
										  .op(1'b0),
										  .A(A_j),
										  .r(r_int),
										  .x_hat(nxt_xhat_j),
										  .nxt_r(r_out),
										  .ready(r_nxt_done));			

/////// Instantiate qadd to calculate dxj = nxt_xhat_j - xhat_j ////////
assign neg_xhat_j[N-1] = ~xhat_j[N-1];
assign neg_xhat_j[N-2:0] = xhat_j[N-2:0];
qadd #(.Q(Q), .N(N)) idxj(.a(nxt_xhat_j),
						  .b(neg_xhat_j),
						  .c(dxj));										  
										  
/////// Instantiate two max modules to compute max_dxj_out and max_xj_out ///////
/*always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		calc_max1 <= 1'b0;
	else
		calc_max1 <= calc_max;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		a_xj <= 'h0;
	else if(start)
		a_xj <= max_xj_in;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		b_xj <= 'h0;
	else if(calc_max1)
		b_xj <= nxt_xhat_j;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		a_dxj <= 'h0;
	else if(start)
		a_dxj <= max_dxj_in;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		b_dxj <= 'h0;
	else if(calc_max1)
		b_dxj <= dxj;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_xj_out <= 'h0;
	else if(update_max)
		max_xj_out <= nxt_max_xj_out;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_dxj_out <= 'h0;
	else if(update_max)
		max_dxj_out <= nxt_max_dxj_out;
end*/

max #(.Q(Q), .N(N)) imax_xj(.a(max_xj_in), .b(nxt_xhat_j), .max(max_xj_out), .which());
max #(.Q(Q), .N(N))	imax_dxj(.a(max_dxj_in), .b(dxj), .max(max_dxj_out), .which());
									  
////// Infer the state flop ///////
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		state <= IDLE;
	else
		state <= nxt_state;
end
		
////// Next state Logic //////
always_comb begin
	
	start_r_int = 1'b0;
	start_r_nxt = 1'b0;
	start_rho = 1'b0;
	start_xhat = 1'b0;
	set_done = 1'b0;
	clr_done = 1'b0;
	//calc_max = 1'b0;
	//update_max = 1'b0;
	rst_init = 1'b1;								// Deasserted by default
	
	case(state)
	
	IDLE:	begin
				rst_init = 1'b0;
				if(start) begin
					clr_done = 1'b1;
					nxt_state = CHECK;
				end
				else
					nxt_state = IDLE;
			end
			
	CHECK: 	begin
				if(A_norm2_j[N-2:0] == 0) begin
					set_done = 1'b1;
					nxt_state = IDLE;
				end
				else begin
					start_r_int = 1'b1;
					nxt_state = CALC_R_INT;
				end
			end
			
	CALC_R_INT: begin
					if(r_int_done) begin
						start_rho = 1'b1;
						nxt_state = CALC_RHO;
					end
					else
						nxt_state = CALC_R_INT;
				end
				
	CALC_RHO:	begin
					if(rho_done) begin
						start_xhat = 1'b1;
						nxt_state = CALC_X;
					end
					else
						nxt_state = CALC_RHO;
				end
			
	CALC_X:		begin
					if(xhat_done) begin
						start_r_nxt = 1'b1;
						//calc_max = 1'b1;
						nxt_state = CALC_R_OUT;
					end
					else 
						nxt_state = CALC_X;
				end
				
	CALC_R_OUT:	begin
					if(r_nxt_done) begin
						//update_max = 1'b1;
						nxt_state = WAIT;
					end
					else
						nxt_state = CALC_R_OUT;
				end
				
	default: 	begin											// This is the WAIT state
					set_done = 1'b1;
					nxt_state = IDLE;
				end
				
	endcase
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		done <= 1'b0;
	else if(clr_done)
		done <= 1'b0;
	else if(set_done)
		done <= 1'b1;
end

endmodule