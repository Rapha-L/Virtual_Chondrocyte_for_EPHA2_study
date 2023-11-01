%% Old verison of script to visualize the result of the timeseries simulation

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
%%

%% Example case
% nodes= [51,54,55,56,58,34,59,23];
% inputvalues = [1,1,1,1,0,1,0,1];
% initial_attr = 1; %sox9+
% inSilicoConditions(2/3,1,'attractorAC_WT.mat',initial_attr,nodes,inputvalues)

%% load case from external simulation
load('inSilico_condition (51 62)_1.mat')

%% find case with transition
for i=1:100
    if condition(1).time_evolution{1}{i}(end,9,1) == 1
        disp(i)
    end
end
%32, 55,66,86

%% Plot timeseries
figure(1)
timeseries = condition(1).time_evolution{1}{32}; % also check 32, 55,66,86
%t=[[1:400],[800:length(timeseries{1})]];
plot([1:length(timeseries)],timeseries(:,[9,10,51,24,20,13],1))
legend('Runx2','Sox9','Cytokines','MMP13','Col2','Col10')

%%
figure(2)
timeseries = condition(1).time_evolution{1}{32};
x = 1:length(timeseries);
timeseries2 = condition(1).time_evolution{2}{38};
x2 = 1:length(timeseries2);
plot(x,timeseries(:,9,1),'b-',x,timeseries(:,10,1),'b--',x,timeseries(:,24,1),'b:',x,timeseries(end,30,1),'b-.',x2,timeseries2(:,9,1),'r-',x2,timeseries2(:,10,1),'r--',x2,timeseries2(:,24,1),'r:',x2,timeseries2(end,30,1),'r-.');
legend('Runx2','Sox9','MMP13','Col2','CRunx2','CSox9','CMMP13','CCol2')

%plot([1:length(timeseries{1})],timeseries{1}(:,[9,10,24,20,13],1))
%legend('Runx2','Sox9','MMP13','Col-II','Col-10')


% figure(2)
% plot([1:length(timeseries{2})],timeseries{2}(:,[9,10,24,20,13],1))
% legend('Runx2','Sox9','MMP13','Col-II','Col-10')