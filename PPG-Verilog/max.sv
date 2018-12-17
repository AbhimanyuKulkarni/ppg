`timescale 1ns / 1ns

module max#(parameter Q = 15,
			parameter N = 32)(
		    input[N-1:0] a,
		    input[N-1:0] b,
		    output reg[N-1:0] max,
			output reg which);
			
always_comb begin

	if((a[N-1] == 1) && (b[N-1] == 0)) begin						// If a is negative and b is positive
		max = b;
		which = 1'b1;
	end
	else if((b[N-1] == 1) && (a[N-1] == 0)) begin					// If a is positive and b is negative
		max = a;
		which = 1'b0;
	end
	else if((a[N-1] == 0) && (b[N-1] == 0)) begin					// If both a and b are positive
		max = (a[N-2:0] > b[N-2:0]) ? a : b;
		which = (a[N-2:0] > b[N-2:0]) ? 1'b0 : 1'b1;
	end
	else begin														// If both a and b are negative
		max = (a[N-2:0] > b[N-2:0]) ? b : a;						// Remember - format is sign and absolute magnitude
		which = (a[N-2:0] > b[N-2:0]) ? 1'b1 : 1'b0;
	end	
end

endmodule

