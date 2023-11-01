%% Write the result of perresult in xl sheet to analyse the percentage of transition per node
% This file is part of the repository Virtual_Chondrocyte_EPHA2_study
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

%% Definition of variables and parameters
n=62; %number of nodes
a=1; % number of the starting attractor for the transition (SOX9=attr2)

%NB: for run2 and run3: change the excel sheet number to write in the excel file 'transition_pernode'

%% loading files
%load the mat containing 'perresult'
load('single_perturbation1sat0.66667.mat') %'trlikrun.....mat'
%load name of nodes (to adapt in the same time as 'n' if necessary)
load('component.mat','componentnames');
names=[componentnames{1,:}];
%name of the excel file? where data are written and saved
xl_file='EphA2_transition_pernode_AC_fromRSox9.xls';

%% wrinting in excel

for c=[1,2]
    if (c==1)
        perturbation_type="d";
    elseif (c==2)
        perturbation_type="u";
    end
    sheet=c; % c for run1 %c+2 for run2  %c+4 for run3
    xlswrite(xl_file,perturbation_type,sheet,'A1') %if run1: c; run2: c+2; run3: c+4
    xlswrite(xl_file,["attr1", "attr2", "attr3"],sheet,'B2') %if run1: c; run2: c+2; run3: c+4
    %write a colonne with the name of node in excel sheet
    xlswrite(xl_file,names',sheet,'A3');%if run1: c; run2: c+2; run3: c+4
    for i=1:n
        j=i+2; % keep few lines for headers above the table
        
        %write value of transition from attractor 'a' to others for each node
        %for perturbation 'c'
        if (~isempty(perresult{a,i,c}))
            xlswrite(xl_file,perresult{a,i,c},sheet,['B' num2str(j)]); %if run1: c; run2: c+2; run3: c+4
        else
             xlswrite(xl_file,[0 0 0],sheet,['B' num2str(j)]);
        end
        
    end
end