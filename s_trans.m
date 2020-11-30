function [s_out]= s_trans(n)
 
  states = (ones(n,1) * 5);
  delays = inf*ones(n,1);
  s_out=[states delays];

