`timescale 1ns / 1ns

module Iterate#(parameter I = 20,
				parameter J = 240,
				parameter Q = 15,
				parameter N = 32)(
				input clk,
				input rst_n,
				input start,
				input[N-1:0] xhat_in[0:J-1],
				input[N-1:0] A[0:I-1][0:J-1],
				input[N-1:0] A_norm2_in[0:J-1],
				input[N-1:0] max_xj_in,
				input[N-1:0] r_in[0:I-1],
				input[N-1:0] lambda,
				output[N-1:0] r_out[0:I-1],
				output[N-1:0] max_xj_out,
				output[N-1:0] max_dxj_out,
				output[N-1:0] xhat_out[0:J-1],
				output reg done);
				
typedef enum logic[1:0] {IDLE, UPDATE, RENEW, CHECK} state_t;
state_t state, nxt_state;

reg update;
reg renew;
reg set_done, clr_done;
reg inc_count;
reg[15:0] count;
reg[N-1:0] max_xj, max_dxj;
reg[N-1:0] r[0:I-1];
reg[N-1:0] xhat[0:J-1];
wire[N-1:0] r_new[0:I-1];
wire[N-1:0] max_xj_new, max_dxj_new;
wire[N-1:0] xhat_new;
wire complete;
wire cnt_complete;

wire[N-1:0] col_A[0:J-1][0:I-1];

genvar c,d;

generate
// Generate col_A, columns of A are rows of col_A
	for(c = 0; c < J; c = c + 1) begin: GEN_COL_A_0
		for(d = 0; d < I; d = d + 1) begin: GEN_COL_A_1
			assign col_A[c][d] = A[d][c];
		end
	end
endgenerate

// Instantiate Update_Element
Update_Element #(.I(I),
				 .Q(Q), 
				 .N(N)) iUpdate(.clk(clk),
								.rst_n(rst_n),
								.start(update), 
								.xhat_j(xhat_in[count]),
								.A_norm2_j(A_norm2_in[count]), 		
								.lambda(lambda),
								.A_j(col_A[count]),
								.max_xj_in(max_xj),
								.max_dxj_in(max_dxj),
								.r_in(r),
								.r_out(r_new),
								.max_xj_out(max_xj_new),
								.max_dxj_out(max_dxj_new),
								.nxt_xhat_j(xhat_new),
								.done(complete));
						
// Counter for counting till J
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		count <= 0;
	else if(start)
		count <= 0;
	else if(inc_count)
		count <= count+1;
end

assign cnt_complete = (count == J);

// Holding registers for max_xj and max_dxj
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_xj <= 0;
	else if(start)
		max_xj <= max_xj_in;
	else if(renew)
		max_xj <= max_xj_new;
end

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		max_dxj <= 0;
	else if(start)
		max_dxj <= 0;
	else if(renew)
		max_dxj <= max_dxj_new;
end

// Holding register for r
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		r <= '{default:'b0};
	else if(start)
		r <= r_in;
	else if(renew)
		r <= r_new;
end

// Holding register for xhat
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		xhat <= '{default:'b0};
	else if(start)
		xhat <= xhat_in;
	else if(renew)
		xhat[count] <= xhat_new;
end

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
	update = 1'b0;
	renew = 1'b0;
	inc_count = 1'b0;
	set_done = 1'b0;
	clr_done = 1'b0;
	
	case(state)
		IDLE:	begin
					if(start) begin
						clr_done = 1'b1;
						update = 1'b1;
						nxt_state = UPDATE;
					end
					else 
						nxt_state = IDLE;
				end
		
		UPDATE:	begin
					if(complete) begin
						renew = 1'b1;
						nxt_state = RENEW;
					end
					else
						nxt_state = UPDATE;
				end
		
		RENEW:	begin
					inc_count = 1'b1;
					nxt_state = CHECK;
				end
		
		CHECK: 	begin
					if(cnt_complete) begin
						set_done = 1'b1;
						nxt_state = IDLE;
					end
					else begin
						update = 1'b1;
						nxt_state = UPDATE;
					end
				end
			
		default: nxt_state = IDLE;
	endcase
end

// S-R Flip Flop for done signal

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		done <= 1'b0;
	else if(clr_done)
		done <= 1'b0;
	else if(set_done)
		done <= 1'b1;
end

assign r_out = r_new;										// All these must be sampled only when done is asserted
assign max_xj_out = max_xj;
assign max_dxj_out = max_dxj;
assign xhat_out = xhat;

endmodule