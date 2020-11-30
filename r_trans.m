function [r_out]= r_trans(n)
 
  states = (ones(n,1) * 6);
  delays = inf*ones(n,1);
  
  r_out=[states delays];

