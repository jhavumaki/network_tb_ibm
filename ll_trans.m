function [ll_out]= ll_trans(lltimes, n, timestep)
  
persistent g;

if isempty(g)
g=0;
end


states = (ones(n,1) * 2);
delay = lltimes((g+1):(g+n));
g=g+n;

ll_out = [states delay];