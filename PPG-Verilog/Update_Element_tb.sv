	`timescale 1ns / 1ns

module Update_Element_tb();

reg clk;
reg rst_n;
reg start;
reg[7:0] xhat_j;
reg[7:0] A_norm2_j;
reg[7:0] lambda;
reg[7:0] A_j[0:9];
reg[7:0] max_xj_in;
reg[7:0] max_dxj_in;
reg[7:0] r_in[0:9];
wire[7:0] r_out[0:9];
wire[7:0] max_dxj_out;
wire[7:0] max_xj_out;
wire[7:0] nxt_xhat_j;
wire done;

// Initialize DUT

Update_Element #(.I(10),
				 .Q(3),
				 .N(8)) iDUT(.clk(clk),
							 .rst_n(rst_n),
							 .start(start),
							 .xhat_j(xhat_j),
							 .A_norm2_j(A_norm2_j),
							 .lambda(lambda),
							 .A_j(A_j),
							 .max_xj_in(max_xj_in),
							 .max_dxj_in(max_dxj_in),
							 .r_in(r_in),
							 .r_out(r_out),
							 .max_dxj_out(max_dxj_out),
							 .max_xj_out(max_xj_out),
							 .nxt_xhat_j(nxt_xhat_j),
							 .done(done));
							 
initial begin
clk = 1'b0;
rst_n = 1'b0;
start = 1'b0;
xhat_j = 8'h00;
A_norm2_j = 8'h70;
lambda = 8'b0_000_0010;
max_xj_in = 8'h00;
max_dxj_in = 8'h00;
A_j = '{8'h24,8'h09,8'h0D,8'h65,8'h01,8'h76,8'hED,8'hF9,8'hC5,8'hE5};
$readmemh("r_initial.hex", r_in);

@(posedge clk);
@(negedge clk);

rst_n = 1'b1;
start =1'b1;

@(posedge clk);
@(negedge clk);

start = 1'b0;

repeat(150) @(posedge clk);
if(done) begin
	$display("YAHOO! Update_Element done!");
	$stop();
end
else begin
	$display("Update_Element not done after 150 cycles.");
	$stop();
end

end

always
#5 clk = ~clk;						
endmodule
						