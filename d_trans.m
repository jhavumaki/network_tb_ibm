function [d_out]= d_trans(life_expectancy,n, timestep)

  death_p = 1.0-exp(-(life_expectancy)/(12/timestep));
  
  delays = (geornd(death_p, n,1) + 1);
  
 d_out = [delays];