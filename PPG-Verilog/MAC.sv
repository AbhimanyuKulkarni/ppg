`timescale 1ns / 1ns

module MAC#(parameter J = 240,
			parameter Q = 15,
			parameter N = 32)(
			input[N-1:0] in_A[0:J-1],
			input[N-1:0] in_B[0:J-1],
			input start,
			input clk,
			input rst_n,
			output[N-1:0] out_C,
			output reg done);
			
typedef enum logic[1:0] {IDLE, MULT, ACCUM} state_t;
state_t state, nxt_state;			
			
wire complete, ovfl;

reg[15:0] i;							// i has to count till J. I think 16 bits should be enough for that
wire i_lt_J;							// Flag to see if i is less than J

// Outputs of state machine
reg set_done, clr_done;
reg mult, accum;
reg clr_i, inc_i;
  
reg[N-1:0] C;  
wire[N-1:0] nxt_C, b_add;
reg[N-1:0] sto_reg;
wire[N-1:0] nxt_sto_reg;

// Instance of the Multiplier
qmults #(.Q(Q), .N(N)) iMULT(.i_multiplicand(in_A[i]), .i_multiplier(in_B[i]), .i_start(mult), .i_clk(clk), .rst_n(rst_n), 
							 .o_result_out(nxt_C), .o_complete(complete), .o_overflow(ovfl));

always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		C <= 'h0;
	else
		C <= nxt_C;
end
							 
// Accumulator storage register
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		sto_reg <= 0;
	else if(start)
		sto_reg <= 0;
	else
		sto_reg <= nxt_sto_reg;
end

assign out_C = sto_reg;
							 
// Instance of the adder								   
qadd #(.Q(Q), .N(N)) iADD(.a(sto_reg), .b(b_add), .c(nxt_sto_reg));
assign b_add = accum ? C : 0;

// Counter for i
always_ff@(posedge clk or negedge rst_n) begin
	if(!rst_n)
		i <= 16'h0000;
	else if(clr_i)
		i <= 16'h0000;
	else if(inc_i)
		i <= i + 1;
end

assign i_lt_J = (i < J-1);													// J-1 because it actually executes one more mult after deasserting

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
				if(i_lt_J) begin
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