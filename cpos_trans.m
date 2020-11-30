function [c_out]= cpos_trans(recoverytimes, deathtimes, txtimes, n, timestep)

  states=[];
  delays=[];

persistent c;

if isempty(c)
c=0;
end

 state_labels = [6 4 3]; %r,d,t
recovery_time = recoverytimes((c+1):(c+n));
death_time= deathtimes((c+1):(c+n));%death times associated with ATB . (deathtimes((c+1):(c+n))).*inf;%deathtimes((c+1):(c+n));
treatment_time = txtimes((c+1):(c+n));
c=c+n;


 
 for i = 1:n
      state = find(df(i,:) == min(df(i,:))); %which indices have min delay time per row
    if (length(state) > 1) % if there are multiple mins, take sample of only 1
        state = randsample(state,1);
    end    
      time = df(i,state); %take actual time
      states(i) = state_labels(state);
      delays(i) = time;
 end   
 
   

  c_out =[states delays];