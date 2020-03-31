`timescale 1ns / 1ps

module testbench();
reg clk;
reg sw;
wire led;

miniCPU cpu(.clk(clk),.sw(sw),.led(led));

initial begin
sw = 0;
clk = 1;
#200 sw = 1;
end

always begin
#5 clk=~clk;
end

endmodule