`timescale 1ns / 1ns

module lasso_tb();

reg clk;
reg rst_n;
reg start_lasso, start_dot;
reg[7:0] y[0:9];
reg[7:0] A[0:9][0:1];
reg[7:0] lambda;
reg[7:0] tol;
wire[7:0] xhat[0:1];
wire done_lasso, done_dot;

reg[7:0] x[0:1];

reg[3:0] c;
reg[7:0] d;
// Initialize DUT

lasso #(.I(10),
		.J(2),
		.Q(1),
		.N(8)) iDUT(.clk(clk),
				   .rst_n(rst_n),
				   .start(start_lasso),
				   .y(y),
				   .A(A),
				   .lambda(lambda),
				   .tol(tol),
				   .xhat(xhat),
				   .done(done_lasso));

// Initialize a dot for calculating y
dot #(.I(10),
	  .J(2),
	  .Q(1),
	  .N(8)) idot(.clk(clk),
				  .rst_n(rst_n),
				  .in_A(A),
				  .in_B(x),
				  .start(start_dot),
				  .out_C(y),
				  .done(done_dot));				   
				   
initial begin
	clk = 1'b0;
	rst_n = 1'b0;
	start_dot = 1'b0;
	start_lasso = 1'b0;
	lambda = 8'b0_000_0001;								// lambda is 0.125
	tol = 8'b0_000_0001;								// tol is 0.0625
	x = '{default: 8'b0_110_0010};						// x = [5.125 5.125 ...]

	for(c = 0; c < 10; c = c + 1) begin
		for(d = 0; d < 2; d = d + 1) begin
			A[c][d] = $random;
		end
	end

	@(posedge clk);
	@(negedge clk);
	
	rst_n = 1'b1;
	start_dot = 1'b1;
	
	@(negedge clk);
	start_dot = 1'b0;
	
	repeat(1000) @(posedge clk);
	
	if(!done_dot)
		$display("Dot product y = Ax not computed after 1000 cycles.");
	else begin
		$display("Dot product computed successfully.");
		
		@(posedge clk);
		@(negedge clk);
		
		start_lasso = 1'b1;
		
		@(negedge clk);
		start_lasso = 1'b0;
		
		repeat(100000) @(posedge clk);
		
		if(!done_lasso)
			$display("Lasso not done after 10k cycles");
		else begin
			$display("Lasso done successfully");
			$display("%p", xhat);
		end
	end
		
	$stop();
		
end

always
#5 clk = ~clk;

endmodule