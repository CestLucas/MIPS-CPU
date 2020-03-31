# MIPS-CPU
An implementation of a MIPS CPU

See Description.docx

In this project, we programmed a five-stage MIPS pipelining CPU using Verilog and downloading it to a ZYBO-7000 development FPGA. 
We use switch 0 on the ZYBO board to start fetching the instruction and use LED 0 to report when processing is done. 
We tested out five R-TYPE instructions to examine our program and added stalls to prevent potential data hazard.
