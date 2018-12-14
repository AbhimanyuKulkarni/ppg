`timescale 1ns / 1ns

/////////////////////////////////////////////////////////
////// This module takes a Ix1 vector and computes //////
////// the sum of squares of all elements. //////////////
/////////////////////////////////////////////////////////

module Square_Sum#(parameter I = 20,
				   parameter Q = 15,
				   parameter N = 32)(
				   input clk,
				   input rst_n,
				   input start,
				   input[N-1:0] in[0:I-1],
				   output[N-1:0] out,
				   output reg done);
			
typedef enum logic[1:0] {IDLE, MULT, ACCUM} state_t;
state_t state, nxt_state;			
			
wire complete, ovfl;

reg[15:0] i;							// i has to count till I. I think 16 bits should be enough for that
wire i_lt_I;							// Flag to see if i is less than I

// Outputs of state machine
reg set_done, clr_done;
reg mult, accum;
reg clr_i, inc_i;
  												
wire[N-1:0] nxt_out, out_add;
reg[N-1:0] to_add;
reg[N-1:0] store_sum;
wire[N-1:0] nxt_sum;

// Instance of the Multiplier
qmults #(.Q(Q), .N(N)) iMULT(.i_multiplicand(in[i]), .i_multiplier(in[i]), .i_start(mult), .i_clk(clk), .rst_n(rst_n), 
							 .o_result_out(nxt_out), .o_complete(complete), .o_overflow(ovfl));


// Sum storage register
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		store_sum <= 0;
	else if(start)
		store_sum <= 0;
	else
		store_sum <= nxt_sum;
end							 

assign out = store_sum;
							 
// Instance of the adder	
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		to_add <= 'h0;
	else
		to_add <= nxt_out;
end
							   
qadd #(.Q(Q), .N(N)) iADD(.a(store_sum), .b(out_add), .c(nxt_sum));
assign out_add = accum ? to_add : 0;

// Counter for i
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		i <= 16'h0000;
	else if(clr_i)
		i <= 16'h0000;
	else if(inc_i)
		i <= i + 1;
end

assign i_lt_I = (i < I-1);													// J-1 because it actually executes one more mult after deasserting

// SR Flip Flop for done
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		done <= 1'b0;
	else if(clr_done)														// Preference to clr_done over set_done
		done <= 1'b0;
	else if(set_done)
		done <= 1'b1;
end

// Infer a state Flop

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		state <= IDLE;
	else
		state <= nxt_state;
end

// Next state and Output logic for FSM

always_comb begin
	
	// Assign default values to outputs 
	mult = 1'b0;	
	accum = 1'b0;
	clr_i = 1'b0;
	inc_i = 1'b0;
	set_done = 1'b0;
	clr_done = 1'b0;
	nxt_state = IDLE;
	
	case(state)
	
	IDLE : 	begin
				if(start) begin
					clr_i = 1'b1;
					clr_done = 1'b1;
					mult = 1'b1;
					nxt_state = MULT;
				end
				else
					nxt_state = IDLE;
			end
			
	MULT : 	begin
				if(complete) begin
					accum = 1'b1;
					nxt_state = ACCUM;
				end
				else 
					nxt_state = MULT;
			end
			
	ACCUM : begin
				if(i_lt_I) begin
					inc_i = 1'b1;
					mult = 1'b1;
					nxt_state = MULT;
				end
				else begin
					set_done = 1'b1;
					nxt_state = IDLE;
				end
			end
			
	default : 	nxt_state = IDLE;
	
	endcase
	
end

endmodule						   

				   