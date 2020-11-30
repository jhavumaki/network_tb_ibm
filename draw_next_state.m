
function out= draw_next_state(curr_state, t, timestep, fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes,deathtimes, monthly_rates, burnin_t, treatment_duration, ipt_duration)  

switch curr_state
case 5
      out = s_trans(1);
case 1
      out = fl_trans(eps, 1, timestep, fltimes, monthly_rates);
case 0
      out = ll_trans(lltimes, 1, timestep);
case 2
      out = cpos_trans(recoverytimes, ATBdeathtimes, txtimes, 1, timestep);
case 3
      out = tx_trans(treatment_duration, 1, timestep);
case 4
      out = [4 inf];
case 6 
    out = r_trans(1);
      
end

