% Script to define conditions and run in silico experiments

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
%%


condition = struct;
newdir= "additionnal simulation";

filename = fullfile(newdir,'inSilico_Epha2_inflam.xls');

%% Define initial state
% From SOX9
attr_init= 1; %sox9+

condition(1).nodes = [62,51,54]; 
condition(1).inputvalues = [0,1,1];
condition(1).initial_attr= attr_init ; %Sox9+

condition(2).nodes = [62,51]; 
condition(2).inputvalues = [0,1];
condition(2).initial_attr= attr_init ; %Sox9+

condition(3).nodes = [62,14]; 
condition(3).inputvalues = [1,1];
condition(3).initial_attr= attr_init ; %Sox9+

condition(4).nodes = [14]; 
condition(4).inputvalues = [1];
condition(4).initial_attr= attr_init ; %Sox9+

condition(5).nodes = [62,14]; 
condition(5).inputvalues = [0,1];
condition(5).initial_attr= attr_init ; %Sox9+

% from Runx2+
attr_init = 2;

condition(6).nodes = [62,51]; 
condition(6).inputvalues = [0,0];
condition(6).initial_attr= attr_init ; %Runx2+

condition(7).nodes = [62,54]; 
condition(7).inputvalues = [0,1];
condition(7).initial_attr= attr_init ; %Runx2+

condition(8).nodes = [54]; 
condition(8).inputvalues = [1];
condition(8).initial_attr= attr_init ; %Runx2+

condition(9).nodes = [62,14]; 
condition(9).inputvalues = [1,1];
condition(9).initial_attr= attr_init ; %Runx2+

condition(10).nodes = [62,14]; 
condition(10).inputvalues = [1,1];
condition(10).initial_attr= attr_init ; %Runx2+

condition(11).nodes = [14]; 
condition(11).inputvalues = [1];
condition(11).initial_attr= attr_init ; %Runx2+

% condition(21).nodes = [62,54,51]; 
% condition(21).inputvalues = [1,1,1];
% condition(21).initial_attr= attr_init ; %Runx2+

for c = 1 : length(condition)
    c
results_runs = Save_insilico_conditions_with_stats(condition(c).nodes,condition(c).inputvalues,condition(c).initial_attr,c,filename);
condition(c).result_transition = results_runs;
end

save(fullfile(newdir,'inSilico_EphA2_Scenarios'),'condition')

function [table2save]= Save_insilico_conditions_with_stats(nodes,inputvalues,initial_attr,sheetNumber,filename)
runs = 3;
results_runs = zeros(runs,2); %nb of attractor
for run = 1 : runs
    run
    [TransitionResult,~] = inSilicoConditions_return(2/3,run,'attractorAC_WT_EPHA2.mat',initial_attr,nodes,inputvalues);
    results_runs(run,:) = TransitionResult ;
end

results_runs(4,:) = mean(results_runs);
results_runs(5,:) = std(results_runs);

colNames = {'Sox9','Runx2'}; % order depends on attractor WT. 
rowNames = {'run1','run2','run3','average','std'};
table2save = array2table(results_runs, 'VariableNames', colNames,'RowNames',rowNames);
writetable(table2save,filename,'WriteVariableNames',true,'WriteRowNames',true,'Sheet',sheetNumber,'Range','A2'); 

xlswrite(filename,[nodes;inputvalues],sheetNumber,'F2');
xlswrite(filename,['init attractor: ' num2str(initial_attr)],sheetNumber,'A9');

end