# network_tb_ibm
Individual based network model of tuberculosis

This model is meant to be run in parallel on a high performance computing cluster. 
The batch script provides the 'ARRAYID' for all files. To pre-generate networks, run “network_generation.R”. 
The main individual-based model file (“main_network_tb_passive.m”) loads a given network depending on the parameter set. 
Interventions are as follows: ‘comm4’ is communitywide-ACF; ‘hh’ is household contact tracing; ‘commCT’ is community contact tracing.
