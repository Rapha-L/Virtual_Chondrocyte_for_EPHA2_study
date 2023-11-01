%% Screen the effect of activation or inhibition of each node/variable of the system, in combinaition with EPHA2 inhibition, on the hypertrophic chondrocyte (Runx2+)
%This function save the simulation result in a .mat file

%% This file is part of the repository Virtual_Chondrocyte_EPHA2_study
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2023 - KU Leuven

%File author(s): R. Lesage

%Distributed under the GPLv3 License.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>

%% Create the structure with input in silico conditions
condition = struct;
n_tot = 61; % All nodes but EPHA2


for c = 1:n_tot  %c = 1:n_tot
 
condition(c).nodes = [62,c]; 
condition(c).inputvalues = [0,1];
condition(c).initial_attr=2 ; 
condition(c).initial_name='Runx2+' ;

condition(c+n_tot).nodes = [62,c]; 
condition(c+n_tot).inputvalues = [0,0];
condition(c+n_tot).initial_attr=2 ;
condition(c+n_tot).initial_name='Runx2+' ;


end



%% Screening simulation
for c = 1 : length(condition)
    c
[results_runs,New_attractor] = Save_insilico_conditions(condition(c).nodes,condition(c).inputvalues,condition(c).initial_attr);
condition(c).result_transition = results_runs;
condition(c).New_attractor = New_attractor;
end

save('inSilico_EphA2-_screen_initRUNX2','condition')

function [results_run,New_attractor]= Save_insilico_conditions(nodes,inputvalues,initial_attr)

    [TransitionResult,newattr] = inSilicoConditions_return(2/3,1,'attractorAC_WT_EPHA2.mat',initial_attr,nodes,inputvalues);
   results_run = TransitionResult;
   New_attractor = newattr;
end