%fast laten transition calcualtion
function [fl_out]= fl_trans(e_rates, n, timestep, fltimes, monthly_rates)
states=[];
persistent k;

if isempty(k) 
k=0;
end

%divide each annual rate (i) by 12 and make array with rate repeating 12 times and append to existing rates
 fast_progression_prob = 1-sum(monthly_rates);%exp(-sum(monthly_rates)); %risk of progression to LL

  fast_prog = double(rand(n,1) > fast_progression_prob); % for each n see if random number drawn from uniform is less than fast prog pro 
   
  states = zeros(1,n);
  times = 5*(12/timestep)*ones(n,1)';
  
  states(fast_prog==1) = 2; %replace  L with active TB is fast progression_prob equals 1
  

  times(fast_prog==1) = fltimes((k+1):(k+sum(fast_prog)));
  k=k+sum(fast_prog);
  
  delay=times;
  fl_out=[states delay];
