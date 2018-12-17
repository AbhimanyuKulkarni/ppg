`timescale 1ns/1ns

module Iterate_tb();

reg clk;
reg rst_n;
reg start;
reg[7:0] xhat_in[0:1];
reg[7:0] A[0:9][0:1];
reg[7:0] A_norm2_in[0:1];
reg[7:0] max_xj_in;
reg[7:0] r_in[0:9];
reg[7:0] lambda;
wire[7:0] r_out[0:9];
wire[7:0] max_xj_out;
wire[7:0] max_dxj_out;
wire[7:0] xhat_out[0:1];
wire done;

Iterate#(.I(10),
		 .J(2),
		 .Q(3),
		 .N(8)) iDUT(.clk(clk),
				.rst_n(rst_n),
				.start(start),
				.xhat_in(xhat_in),
				.A(A),
				.A_norm2_in(A_norm2_in),
				.max_xj_in(max_xj_in),
				.r_in(r_in),
				.lambda(lambda),
				.r_out(r_out),
				.max_xj_out(max_xj_out),
				.max_dxj_out(max_dxj_out),
				.xhat_out(xhat_out),
				.done(done));

initial begin

	clk = 1'b0;
	rst_n = 1'b0;
	start = 1'b0;
	xhat_in = '{default:8'h00};
	A[0] = '{8'h24,8'h81};
	A[1] = '{8'h09,8'h63};
	A[2] = '{8'h0D,8'h8D};
	A[3] = '{8'h65,8'h12};
	A[4] = '{8'h01,8'h0D};
	A[5] = '{8'h76,8'h3D};
	A[6] = '{8'hED,8'h8C};
	A[7] = '{8'hF9,8'hC6};
	A[8] = '{8'hC5,8'hAA};
	A[9] = '{8'hE5,8'h77};
	$display(A);
	A_norm2_in = '{8'h22,8'h00};
	max_xj_in = 8'h00;
	r_in = '{8'h44,8'h50,8'h74,8'h0C,8'h08,8'h26,8'hD4,8'hC0,8'hEC,8'h8C};
	lambda = 8'b0_000_0010;											// lambda is 0.125
	
	@(posedge clk);
	@(negedge clk);
	
	rst_n = 1'b1;
	start = 1'b1;
	
	@(posedge clk);
	@(negedge clk);
	
	start = 1'b0;

	repeat(150) @(posedge clk);
	if(!done)
			$display("Iterate not done after 150 cycles");
	else
			$display("YAHOO! Iterate done successfully");
			
	$stop();
	
end

always
#5 clk = ~clk;

endmodule