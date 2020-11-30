function screen_out=community_contact_screening(x, state_df, t, CommConnects)

if ~isempty(CommConnects{x})

    if length(CommConnects{x})> 3
        n=4;
    elseif length(CommConnects{x})< 4
            n=length(CommConnects{x});

    end

commIDs= datasample(CommConnects{x},n,'Replace', false);

commIDs_states = reshape(state_df(commIDs,2),size(commIDs)); %states for all HH members
%S_states =find(commIDs_states==5);
%state_df(commIDs(S_states),2) =  7; %current state is 
%state_df(commIDs(S_states),3) =  5; %next state is 5
%state_df(commIDs(S_states),4) = t+6;
      
L_states =find(commIDs_states==1 | commIDs_states==0);
state_df(commIDs(L_states),2) =  8; %current state is ptL
state_df(commIDs(L_states),3) =  6; %next state is 6 i.e., R
state_df(commIDs(L_states),4) = t+6;
      
C_states =find(commIDs_states==2);
state_df(commIDs(C_states),2) =  3; %current state is ptL
state_df(commIDs(C_states),3) =  6; %next state is 6 i.e., R
state_df(commIDs(C_states),4) = t+6;

latent=commIDs_states==1 | commIDs_states==0;
ATB=commIDs_states==2;
screen_out= [commIDs state_df(commIDs(:),2) state_df(commIDs(:),3) state_df(commIDs(:),4) latent ATB];

elseif isempty(CommConnects{x})
    screen_out=[];
end

end
