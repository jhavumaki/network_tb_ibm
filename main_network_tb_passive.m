%adapted from Kasaie et al. 2014, main passive surveillance only model
function main_network_tb_passive(ARRAYID)

pp=1;
HH = 5; %size of households
HH_contacts = HH-1; %number of contacts per household
init_N= 100000; %total population
maxID = init_N;
IDs = [1:maxID]'; 
numSamples =30000000;% number transition times to pre-draw
load lhs_9000_200110.mat %load natural history parameter
timestep=1; %monthly timesteps, this could be changed to annual i.e., 12
burnin_t = 6600; %time to steady state
tspan = 6672; %total simulation time
totalrandseed=5; %number of stochastic runs per parameter set
 
%range of parameter sets for each cluster batch job (i.e., 10, but can change)
minLHS=10*(ARRAYID-1)+1;
maxLHS=ARRAYID*10;
numLHS=maxLHS-minLHS+1;
 
rdmseeds=csvread('rdmseeds_independant9000by5by10.csv',1); %pre-drawn random seeds to account for stochasticity
 

%initialize record keeping structures
incidence={};

AnnualCInf={};
AnnualHHInf={};
newTx={};



for l = minLHS:maxLHS %range of parameter sets
 
connects={};    
CommConnects={};
 
network=paramLHS(25,l);
array=paramLHS(26,l);

%load network file
networkfile = strcat('Network', num2str(network),'array',num2str(array), '.csv')   
NetArray=csvread(networkfile,1,1); 
connMat=[NetArray(:,1) NetArray(:,2) NetArray(:,3)];
   
%add in missing households that are not connected to other households
    nrConnectedID= length(unique(connMat(:,1)));
    if nrConnectedID(1)<init_N
    missingIDs=setdiff(IDs,unique(connMat(:,1)));
    missingIDs= [missingIDs zeros(length(missingIDs),1) zeros(length(missingIDs),1)];
    connMat=sortrows([connMat;missingIDs]);
    end
    
    cutPt=find(diff(connMat(:,1))>0);
    a=size(connMat);
    connects=mat2cell(connMat, diff([0;(cutPt);a(1)]),3);
 
    dimConn=size(connects);
   
    for i = 1:maxID
    CommConnects{i}=connects{i}(connects{i}(:,3)~=0,2);
    end
    
   
 
%load natural history parameters
      theta = 1/paramLHS(1,l);
      epses = paramLHS(2:6,l);
      tau = paramLHS(7,l);
      gamma =paramLHS(8,l);
      kappa =paramLHS(9,l);
      p_ipt =paramLHS(10,l);
      ATB_treat = paramLHS(11,l); 
      treatment_duration = 12/paramLHS(12,l);
      ipt_duration = 12/paramLHS(13,l);
      betaHH = paramLHS(14,l);
   
      betaComm=paramLHS(29,l);
      
%set to 0 by default, but here is where to switch on imported TB cases
      imported=0;%((paramLHS(14,l)*paramLHS(15,l))/init_N)*paramLHS(24,l);
 
      infection_protection = paramLHS(16:22,l);
    
    t=[];
    
for ss = 1:totalrandseed %

 %as it is none take first 5
rng(rdmseeds(l, ss)); %need to set seed before random number draw
s=rng
    
 
monthly_rates=repelem(epses/(12/timestep), (12/timestep));
 
 
%random number draws of transition times
fltimes = datasample(1:((12/timestep)*5), numSamples,'weights', monthly_rates);
lltimes = geornd(1.0-exp(-tau/(12/timestep)),numSamples,1)+1;
recoverytimes = geornd(1.0-exp(-gamma/(12/timestep)), numSamples,1) + 1;
ATBdeathtimes = geornd(1.0-exp(-kappa/(12/timestep)), numSamples,1) + 1;
txtimes =geornd(1.0-exp(-ATB_treat/(12/timestep)),numSamples,1) + 1;%geornd(1.0-exp(-paramValues(10,2)/(12/timestep)),numSamples,1) + 1;
deathtimes = nan+zeros(numSamples,1);% geornd(1.0-exp(-(paramValues(1,1))/(12/timestep)), numSamples,1) + 1;%
 
%%immunity
inf_multiplier = 1.0 - infection_protection; 
inf_multiplierS =  struct('S', inf_multiplier(1), 'E',inf_multiplier(2),'L',inf_multiplier(3),'C',inf_multiplier(4),'T',inf_multiplier(5),'D',inf_multiplier(6),'R',inf_multiplier(7));
 
%Initial States
numATB=1; %infecitous seed
seed= repmat(5,maxID,1);
seed(randsample(1:maxID,numATB))=2;
 
S_count=sum(seed==5); %susceptible
E_count=sum(seed==1); %early latent
L_count=sum(seed==0); %late latent
C_count=sum(seed==2); %active TB
T_count=sum(seed==3); %treatment
D_count=sum(seed==4); %death
R_count=sum(seed==6); %recovered
PTs_count=sum(seed==7); %preventive therapy
PTl_count=sum(seed==8); %preventive therapy
Init_S=struct('S',S_count, 'E',E_count, 'L',L_count, 'C',C_count, 'T', T_count, 'D',D_count, 'R',R_count, 'PTs', PTs_count, 'PTl', PTl_count);
 
%% Create a cell array to store agent information

%draw next states, times 
init_states=arrayfun(@(x) draw_next_state(x,t, timestep,fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes, deathtimes, monthly_rates, burnin_t, treatment_duration, ipt_duration), seed, 'UniformOutput', false);
 
 
init_states = cat(1,init_states{:});%combine init_states from multiple structs to a single struct
next_state = init_states(:,1);%[init_states.states]';
next_event_t =init_states(:,2);%[init_states.delay]';

% Allocate individuals into HHs 
num_hh = round(maxID/HH) ;
hh_id = ceil(IDs/5);%randsample(1:num_hh,init_N, true)';
 

%to keep track of total number of people per state by time step
h=fieldnames(Init_S)';%headings
M = zeros(tspan+1, length(fieldnames(Init_S))); 
state_ts = [h; num2cell(M)];
state_ts(2,:) = struct2cell(Init_S);
 
death_time= d_trans(theta,init_N, timestep);% sum([Init_S.counts{:}]));
 
%individual tracking 
state_df = [IDs seed next_state next_event_t death_time hh_id];


%initialize record keeping
curr_state =[];
init_states =[];
next_state =[];
next_event_t =[];
hh_id =[];
state_dfo =[];
incATB=[];
atb_pos=[];
cumI=0;
y=1;
ATB_ID=[];
HH_Comm_ATB=[];

%start simulation
for t = 1:tspan

infectedsC=state_df(state_df(:,2)==2); %new infected individuals
ind_foi=zeros(maxID,1); 

%force of infection
if length(infectedsC)>0

%household
hh=unique(ceil([infectedsC]/5)); 
hh_ids=[(hh*5-4)+0 (hh*5-4)+1 (hh*5-4)+2 (hh*5-4)+3 (hh*5-4)+4];
hh_states = reshape(state_df(hh_ids,2),size(hh_ids));
num_C_HH=sum(hh_states == 2,2);
foi_HH = (betaHH * num_C_HH); %num_B_HH +ci_rho*
foinx=size(hh_states);
foivalues = (hh_states~=2  & hh_states~=3 & hh_states~=4 ).*repmat(foi_HH, [1 foinx(2)]);
indices=hh_ids(hh_states~=2  & hh_states~=3 & hh_states~=4);
foi_hh=zeros(init_N,1);
foi_hh(indices)=foivalues(foivalues>0); 
 
%community
allC=cat(1,CommConnects{infectedsC}); %all community people connected to a C; new
 
numC_comm=zeros(init_N,1);
if length(allC)==1
    numC_comm(allC) = 1; 
else 
numC_comm(unique(allC))=hist(allC, unique(allC));
end
    
foi_Comm = (betaComm * (numC_comm));  %numB_comm + ci_rho*
foi= foi_hh+foi_Comm + ones(init_N,1)*imported;
 
statearray=state_df(:,2);
S_class=(statearray==5);
E_class=(statearray==1);
L_class=(statearray==0);
R_class=(statearray==6);
        ind_foi(S_class)  = (1.0-exp(-1*foi(S_class)*inf_multiplierS.S));
        ind_foi(E_class)  = (1.0-exp(-1*foi(E_class)*inf_multiplierS.E));
        ind_foi(L_class)  = (1.0-exp(-1*foi(L_class)*inf_multiplierS.L));     
        ind_foi(R_class)  = (1.0-exp(-1*foi(R_class)*inf_multiplierS.R));     
end 
 
 
%who becomes infected
ind_inf = rand(length(ind_foi),1) < ind_foi;
newLatent=find(ind_inf==1); %new latent infections at this time step
 

cumI=sum(ind_inf)+cumI;
 
% Change next state for newly infected individuals and set transition time to this step
  state_df(ind_inf,3)= 1;
  state_df(ind_inf,4) = t;
  
% Get all individuals who have transitioned on this time step
  transition_ids = state_df(:,4) == t;
  num_transitions = sum(transition_ids);
 
incATB = sum(state_df(transition_ids,3) ==2);%new Active TB
newTx{l-(minLHS)+1, ss}(t)= sum(state_df(transition_ids,3) ==3); %new treatment
 
  if (num_transitions > 0) 
% draw a set of next states for individuals that transitions
      state_df(transition_ids,2) =  state_df(transition_ids,3);
     sum(state_df(transition_ids,2) ==2);
 
      states1=state_df(transition_ids,2);
      curr_time=repmat(t, length(states1),1);
      next_states = arrayfun(@(x) draw_next_state(x,t, timestep,fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes, deathtimes, monthly_rates, burnin_t, treatment_duration, ipt_duration), states1, 'UniformOutput', false);
      next_states = cat(1,next_states{:});
      state_df(transition_ids,3) = next_states(:,1);
      state_df(transition_ids,4) = t + next_states(:,2);
  end
 


  if t>burnin_t || (t==1 && burnin_t==0)
     ATB_ID_TEMP=[]; 
      IDsFound=find(transition_ids==1); %which IDs transitioned
      TxIDs=IDsFound(state_df(IDsFound,2)==3); %which IDs transitioned to Tx
      
%community/HH contact prevalence among detected
state_dfTemp=[state_df;
100001 0 0 Inf 0 20001; 
100002 1 0 Inf 0 20001; 
100003 2 0 Inf 0 20001; 
100004 3 0 Inf 0 20001; 
100005 4 0 Inf 0 20001; 
100006 5 0 Inf 0 20001;
100007 6 0 Inf 0 20001;
100008 7 0 Inf 0 20001;
100009 8 0 Inf 0 20001;]; %adding this in to ensure all states are included in the cross tab below

r=crosstab(state_dfTemp(:,2), state_dfTemp(:,6));
r=r(:, 1:(end-1));
HH_Comm_ATB(t-burnin_t,1)= mean(sum(r(1:2,:),1)/5); %hh contacts latent
HH_Comm_ATB(t-burnin_t,2)= mean(sum(r(3,:),1)/5); %hh contacts active

a=hist(NetArray(:,1), unique(NetArray(:,1)))'; %lengths of community contacts
NetArrayTemp=[NetArray; state_dfTemp(100001:end, 1) state_dfTemp(100001:end, 1) repmat(0,9,1)];

NetArrayTemp=[NetArrayTemp state_dfTemp(NetArrayTemp(:,2),2)];
r=crosstab(NetArrayTemp(:,1),NetArrayTemp(:,4)); 
r=r(1:(end-9),:);
% sum by row for latent (column 1 and 2 correspond to state 0 and 1 and we divide by number in each community contact group and take the mean
HH_Comm_ATB(t-burnin_t,3)= mean(sum(r(:,1:2),2)./a); %community contacts latent
HH_Comm_ATB(t-burnin_t,4)= mean(sum(r(:,3),2)./a); %community contacts active

%keep track of who moved to I and who moved out to calculate infectious period
ATB_IDs=IDsFound(state_df(IDsFound,2)==2);
ATB_ID_TEMP=[ATB_IDs t*ones(length(ATB_IDs),1) state_df(ATB_IDs,3) state_df(ATB_IDs,4)] ;
ATB_ID=[ATB_ID; ATB_ID_TEMP];

 end
  
 %births and deaths 
 disease_death_ids = state_df(:,2) == 4;
 natural_death_ids = state_df(:,5) == t;
 
 if (sum(natural_death_ids) > 0) 
   state_df(natural_death_ids,2) = 5;
   state_df(natural_death_ids,3) = 5;
   state_df(natural_death_ids,4) = inf;
   state_df(natural_death_ids,5) = t+d_trans(theta,sum(natural_death_ids), timestep);
 end
 
%births reshuffling scheme
if sum(disease_death_ids)>0
state_dfsub=state_df(disease_death_ids,1);
for i = 1:length(state_dfsub)
x=[];
 deg1=[];
 deg2=[];
 deg3=[];
id=state_dfsub(i);
x=CommConnects{id};
if isempty(x)
    state_df(id,2)=5;
    state_df(id,3)=5;
    state_df(id,4)=inf;
    state_df(id,5)=t+d_trans(theta,1, timestep);
end
 
if ~isempty(x)
    x=x(state_df(x,2)~=4); %remove contacts who are dead
    x=x(state_df(x,2)~=2); %remove contacts who are currently infected
    if isempty(x)
    state_df(id,2)=5;
    state_df(id,3)=5;
    state_df(id,4)=inf;
    state_df(id,5)=t+d_trans(theta,1, timestep);
    end
    if ~isempty(x)
    deg1=randsample(x,1);
    state_df(id,2:5) =state_df(deg1,2:5);
    end
end
 
 
if isempty(deg1)
        state_df(deg1,2)=5;
        state_df(deg1,3)=5;
        state_df(deg1,4)=inf;
        state_df(deg1,5)=t+d_trans(theta,1, timestep);
end
 
if ~isempty(deg1)
deg1id=deg1;
deg1=CommConnects{deg1};
    deg1=deg1(state_df(deg1,2)~=4); %remove contacts who are dead
    deg1=deg1(state_df(deg1,2)~=2); %remove contacts who are currently infected
    if isempty(deg1)
    state_df(deg1,2)=5;
    state_df(deg1,3)=5;
    state_df(deg1,4)=inf;
    state_df(deg1,5)=t+d_trans(theta,1, timestep);
    end
    if ~isempty(deg1)
        deg2=randsample(deg1,1);
        state_df(deg1id,2:5)=state_df(deg2,2:5); 
    end
end
 
    if isempty(deg2)
            state_df(deg2,2)=5;
            state_df(deg2,3)=5;
            state_df(deg2,4)=inf;
            state_df(deg2,5)=t+d_trans(theta,1, timestep);
    end
    if ~isempty(deg2)
    deg2id=deg2;
    deg2=CommConnects{deg2};
    deg2=deg2(state_df(deg2,2)~=4); %remove contacts who are dead
    deg2=deg2(state_df(deg2,2)~=2); %remove contacts who are currently infected
    if isempty(deg2)   
            state_df(deg2,2)=5;
            state_df(deg2,3)=5;
            state_df(deg2,4)=inf;
            state_df(deg2,5)=t+d_trans(theta,1, timestep);
    end
    if ~isempty(deg2)
            deg3=randsample(deg2,1);
            state_df(deg2id,2:5)=state_df(deg3,2:5); 
            state_df(deg3,2)=5;
            state_df(deg3,3)=5;
            state_df(deg3,4)=inf;
            state_df(deg3,5)=t+d_trans(theta,1, timestep);
      
            end
    end
end
end
 
%update total number of people per state
 state_ts{t+2,1} =sum(state_df(:,2)==5);
 state_ts{t+2,2} =sum(state_df(:,2)==1);
 state_ts{t+2,3} =sum(state_df(:,2)==0);
 state_ts{t+2,4} =sum(state_df(:,2)==2);
 state_ts{t+2,5} =sum(state_df(:,2)==3);
 state_ts{t+2,6} =sum(state_df(:,2)==4);
 state_ts{t+2,7} =sum(state_df(:,2)==6);
 state_ts{t+2,8} =sum(state_df(:,2)==7);
 state_ts{t+2,9} =sum(state_df(:,2)==8);
 
%reset index
if mod(t,12/timestep)==0
m = 12/timestep;
else m = mod(t,12/timestep);
end

% to keep track of community vs. household attributable infections
foi_CommARTI=foi_Comm(newLatent);
foi_HHARTI=foi_hh(newLatent);
 
commOnlyInf = sum(foi_HHARTI==0 & foi_CommARTI>0);
hhOnlyInf = sum(foi_HHARTI>0 & foi_CommARTI==0);
 
foiindx = foi_HHARTI>0 & foi_CommARTI>0;
competingFoi = [foi_CommARTI(foiindx) foi_HHARTI(foiindx)];
if sum(size(competingFoi))>0
commPr = competingFoi(:,1)./(competingFoi(:,1)+competingFoi(:,2));
else commPr=0;
end
 
commInf(m)= sum(commPr< rand(length(commPr),1))+commOnlyInf;
hhInf(m)= sum(commPr> rand(length(commPr),1))+hhOnlyInf;
 
 
deaths(m) = sum(natural_death_ids)+ sum(disease_death_ids);
atb_pos(m) = incATB;
 
 if t/(12/timestep)==round(t/(12/timestep))
incidence{l-(minLHS)+1, ss}(y)=sum(atb_pos); %annual incidence

    
     AnnualCInf{l-(minLHS)+1, ss}(y)= sum(commInf);
     AnnualHHInf{l-(minLHS)+1, ss}(y)=sum(hhInf);
     
     y=y+1; 
    
     tx=state_ts{t+2,5};
     atb_pos=[];
     commInf=[];
     hhInf=[];
     deaths=[];
     incATB=[];
 end

a=[];
infTemp=[];
 
end
end

%go back to beginning of pre-drawn transition times vector

clear fl_transvNoParams2_matPredraw6
clear ll_transvNoParams2_matPredraw 
clear cpos_transvNoParams2_matPredraw4
clear tx_transvNo6
clear d_transvNoParams2_matPredraw2
clear allC
clear numC_comm
 
%save state_ts: number of people per state per time
filename_ts= strcat('state_ts','net', num2str(network),'array',num2str(array), 'lhs', num2str(l),'rng',num2str(rdmseeds(l,ss)),'seed', num2str(ss),'none','.csv')
fid = fopen(filename_ts, 'w'); 
fprintf(fid, '%s,', state_ts{1,:}) ;
fclose(fid);
dlmwrite(filename_ts, cell2mat(state_ts(2:end,:)),'roffset',1, '-append')

%%save prop ATB/LTBI by household and community contacts
filename_ts= strcat('ATB_comm_hh','net', num2str(network),'array',num2str(array), 'lhs', num2str(l),'rng',num2str(rdmseeds(l,ss)),'seed', num2str(ss),'none','.csv');
fid = fopen(filename_ts, 'w') ;
fclose(fid);
dlmwrite(filename_ts, HH_Comm_ATB,'roffset',1, '-append')

%save length of ATB time; infectious period
filename_ts= strcat('ATB_infPer','net', num2str(network),'array',num2str(array), 'lhs', num2str(l),'rng',num2str(rdmseeds(l,ss)),'seed', num2str(ss),'none','.csv');
fid = fopen(filename_ts, 'w') ;
fclose(fid);
dlmwrite(filename_ts, ATB_ID,'roffset',1, '-append')

pp=pp+1;

end
 
end

%save other files 1 per cluster job 
 
 iMatNames=cell(1,numLHS*totalrandseed);
 incidenceMat=[];
 newTxMat=[];


 infPerMat=[];
 commInfMat=[];
 hhInfMat=[];
 k=1;
 kk=1;
 for ss= 1:totalrandseed
         for l = minLHS:maxLHS
    incidenceMat(:,k) = incidence{l-(minLHS)+1 ,ss}';
    newTxMat(:,k) = newTx{l-(minLHS)+1 ,ss}';
    commInfMat(:,k) = AnnualCInf{l-(minLHS)+1, ss}';
    hhInfMat(:,k) =AnnualHHInf{l-(minLHS)+1, ss}';
    iMatNames{k} = char(strcat('lhs',num2str((l)),'seed', num2str(rdmseeds(l,ss)), 'none'));

    k=k+1;
    kk=kk+1;
         end
         
 end
 
 
%annual incidence
  header= iMatNames;
 filename_dfo= strcat('inc_', 'lhs',num2str(minLHS), 'none',date, '.csv'); 
 fid = fopen(filename_dfo, 'w') ;
 fprintf(fid, '%s,', header{:}) ;
 fclose(fid);
 dlmwrite(filename_dfo, incidenceMat,'roffset',1, '-append')
 
 
%community attributable infections
  header= iMatNames;
 filename_dfo= strcat('annualC_', 'lhs',num2str(minLHS), 'none',date, '.csv'); 
 fid = fopen(filename_dfo, 'w') ;
 fprintf(fid, '%s,', header{:}) ;
 fclose(fid);
 dlmwrite(filename_dfo, commInfMat,'roffset',1, '-append')
 

%househod attributable infections
 header= iMatNames;
 filename_dfo= strcat('annualHH_', 'lhs',num2str(minLHS), 'none',date, '.csv'); 
 fid = fopen(filename_dfo, 'w') ;
 fprintf(fid, '%s,', header{:}) ;
 fclose(fid);
 dlmwrite(filename_dfo, hhInfMat,'roffset',1, '-append')
 
 
%new treatment

 header= iMatNames;
 filename_dfo= strcat('newTx_', 'lhs',num2str(minLHS), 'none',date, '.csv'); 
 fid = fopen(filename_dfo, 'w') ;
 fprintf(fid, '%s,', header{:}) ;
 fclose(fid);
 dlmwrite(filename_dfo, newTxMat,'roffset',1, '-append')
 

end

