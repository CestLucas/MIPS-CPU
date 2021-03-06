# 
# Synthesis run script generated by Vivado
# 

set_param xicom.use_bs_reader 1
set_msg_config -id {HDL 9-1061} -limit 100000
set_msg_config -id {HDL 9-1654} -limit 100000
set_msg_config -id {Synth 8-256} -limit 10000
set_msg_config -id {Synth 8-638} -limit 10000
create_project -in_memory -part xc7z010iclg400-1L

set_param project.singleFileAddWarning.threshold 0
set_param project.compositeFile.enableAutoGeneration 0
set_param synth.vivado.isSynthRun true
set_property webtalk.parent_dir {C:/Users/Lucas/Documents/Final Project/Final Project.cache/wt} [current_project]
set_property parent.project_path {C:/Users/Lucas/Documents/Final Project/Final Project.xpr} [current_project]
set_property default_lib xil_defaultlib [current_project]
set_property target_language Verilog [current_project]
set_property ip_output_repo {c:/Users/Lucas/Documents/Final Project/Final Project.cache/ip} [current_project]
set_property ip_cache_permissions {read write} [current_project]
read_verilog -library xil_defaultlib {{C:/Users/Lucas/Documents/Final Project/Final Project.srcs/sources_1/new/Final Project.v}}
foreach dcp [get_files -quiet -all *.dcp] {
  set_property used_in_implementation false $dcp
}
read_xdc {{C:/Users/Lucas/Documents/Final Project/Final Project.srcs/constrs_1/imports/VivadoDemo/ZYBO_MASTER.xdc}}
set_property used_in_implementation false [get_files {{C:/Users/Lucas/Documents/Final Project/Final Project.srcs/constrs_1/imports/VivadoDemo/ZYBO_MASTER.xdc}}]


synth_design -top miniCPU -part xc7z010iclg400-1L


write_checkpoint -force -noxdef miniCPU.dcp

catch { report_utilization -file miniCPU_utilization_synth.rpt -pb miniCPU_utilization_synth.pb }
