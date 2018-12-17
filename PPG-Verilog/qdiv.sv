`timescale 1ns / 1ns
//////////////////////////////////////////////////////////////////////////////////
// Company: 			Burke
// Engineer: 			Tom Burke
// 
// Create Date:		19:39:14 08/24/2011 
// Design Name: 	
// Module Name:		qdiv.v
// Project Name:		Fixed-point Math Library (Verilog)
// Target Devices: 
// Tool versions:		Xilinx ISE WebPack v14.7
// Description: 		Fixed-point division in (Q,N) format
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
//	Revision 0.02 - 25 May 2014
//							Updated to fix an error
//
// Additional Comments: Based on my description on youtube:
//			http://youtu.be/TEnaPMYiuR8
//
//////////////////////////////////////////////////////////////////////////////////
 
module qdiv #(
	//Parameterized values
	parameter Q = 15,
	parameter N = 32
	)
	(
	input 	[N-1:0] i_dividend,
	input 	[N-1:0] i_divisor,
	input 	i_start,
	input 	i_clk,
	input	rst_n, 
	output 	[N-1:0] o_quotient_out,
	output 	o_complete,
	output	o_overflow,
	output reg 	div_0
	);
 
	reg [2*N+Q-3:0]	reg_working_quotient;	//	Our working copy of the quotient
	reg [N-1:0] 		reg_quotient;				//	Final quotient
	reg [N-2+Q:0] 		reg_working_dividend;	//	Working copy of the dividend
	reg [2*N+Q-3:0]	reg_working_divisor;		// Working copy of the divisor
 
	reg [N-1:0] 			reg_count; 		//	This is obviously a lot bigger than it needs to be, as we only need 
													//		count to N-1+Q but, computing that number of bits requires a 
													//		logarithm (base 2), and I don't know how to do that in a 
													//		way that will work for everyone
										 
	reg					reg_done;			//	Computation completed flag
	reg					reg_sign;			//	The quotient's sign bit
	reg					reg_overflow;		//	Overflow flag
	reg 				set_done, clr_done;
	reg 				ongoing;
 
	assign o_quotient_out[N-2:0] = reg_quotient[N-2:0];	//	The division results
	assign o_quotient_out[N-1] = reg_sign;						//	The sign of the quotient
	assign o_complete = reg_done;
	assign o_overflow = reg_overflow;
 
	always @( posedge i_clk or negedge rst_n) begin
		if(!rst_n) begin
			ongoing <= 1'b0;				//	Initial state is to not be doing anything
			reg_overflow <= 1'b0;			//		And there should be no overflow present
			reg_sign <= 1'b0;				//		And the sign should be positive
			reg_working_quotient <= 0;	
			reg_quotient <= 0;				
			reg_working_dividend <= 0;	
			reg_working_divisor <= 0;		
			reg_count <= 0;
			div_0 <= 0;
			clr_done <= 0;
			set_done <= 0;
		end			
		else if( i_start ) begin										//	This is our startup condition
			clr_done <= 1'b1;												//	We're not done			
			div_0 <= 1'b0;
			
			if(i_divisor[N-2:0] == 0) begin
				div_0 <= 1'b1;
				set_done <= 1'b1;
			end
			
			ongoing <= 1'b1;
			reg_count <= N+Q-1;											//	Set the count
			reg_working_quotient <= 0;									//	Clear out the quotient register
			reg_working_dividend <= 0;									//	Clear out the dividend register 
			reg_working_divisor <= 0;									//	Clear out the divisor register 
			reg_overflow <= 1'b0;										//	Clear the overflow register

			reg_working_dividend[N+Q-2:Q] <= i_dividend[N-2:0];				//	Left-align the dividend in its working register
			reg_working_divisor[2*N+Q-3:N+Q-1] <= i_divisor[N-2:0];		//	Left-align the divisor into its working register

			reg_sign <= i_dividend[N-1] ^ i_divisor[N-1];		//	Set the sign bit
			end 
		else if(ongoing) begin
			reg_working_divisor <= reg_working_divisor >> 1;	//	Right shift the divisor (that is, divide it by two - aka reduce the divisor)
			reg_count <= reg_count - 1;								//	Decrement the count

			//	If the dividend is greater than the divisor
			if(reg_working_dividend >= reg_working_divisor) begin
				reg_working_quotient[reg_count] <= 1'b1;										//	Set the quotient bit
				reg_working_dividend <= reg_working_dividend - reg_working_divisor;	//		and subtract the divisor from the dividend
				end
 
			//stop condition
			if(reg_count == 0) begin
				set_done <= 1'b1;										//	If we're done, it's time to tell the calling process
				ongoing <= 1'b0;
				reg_quotient <= reg_working_quotient;			//	Move in our working copy to the outside world
				if (reg_working_quotient[2*N+Q-3:N]>0)
					reg_overflow <= 1'b1;
					end
			else
				reg_count <= reg_count - 1;	
			end
		end
		
	always_ff@(posedge i_clk or negedge rst_n) begin
		if(!rst_n)
			reg_done <= 1'b0;
		else if(clr_done)														// Preference to clr_done over set_done
			reg_done <= 1'b0;
		else if(set_done)
			reg_done <= 1'b1;
	end
		
endmodule