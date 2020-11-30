function [tx_out] = tx_trans(treatment_duration,n, timestep)

  states = (ones(n,1) * 6);
  delays = 6;%(treatment_duration*12/timestep)*ones(n,1);
   

tx_out = [states delays];