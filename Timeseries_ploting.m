%% Updated Script for visualization of timeseries simulation results

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

%% Example case
% nodes= [51,54,55,56,58,34,59,23];
% inputvalues = [1,1,1,1,0,1,0,1];
% initial_attr = 1; %sox9+
% inSilicoConditions(2/3,1,'attractorAC_WT.mat',initial_attr,nodes,inputvalues)

%% load case from external simulation
load('inSilico_condition (62  51)_2.mat')
load('component.mat') %index to name conversion

%% find case with transition
series=timeseries;
for i=1:100
    if series{i}(end,9,1) == 1
        disp(i)
    end
end
%32, 55,66,86

%% plot timeseries
figure(1)
s = 66; %NB: change it a t line 85 too
series = timeseries{s}; % also check 32, 55,66,86
%t=[[1:400],[800:length(timeseries{1})]];
plot([1:length(series)],series(:,[9,10,24,20,13,7],1))

legend('RUNX2','SOX9','MMP13','COL-II','COL-X','B-catenin')
xlabel('Timesteps')
ylabel('Activity profile')
title({'Effect of perturbation activity profiles', 'of selected entities over time'})

Sizes = cellfun('size',timeseries,1);
% max(Sizes)
% min(Sizes)

x_max= min(Sizes);
storeMatrice = zeros(x_max,62,3,100);

for s =1:100
    storeMatrice(:,:,:,s)=timeseries{s}(1:x_max,:,:);
end
%%
GlobalOnlyMatrice=zeros(x_max,62,100);
GlobalOnlyMatrice(:,:,:)= storeMatrice(:,:,1,:);

averages = mean(GlobalOnlyMatrice,3);
stds = std(GlobalOnlyMatrice,0,3);
Plusstd = averages + stds;
Lessstd = averages - stds;





%%
figure(2)
x = 1:x_max;
Nodes=[7,29,9];
plot(x,averages(:,Nodes(1)),'b-',x,Plusstd(:,Nodes(1)),'b--',x,Lessstd(:,Nodes(1)), 'b--', ...
    x,averages(:,Nodes(2)),'r-',x,Plusstd(:,Nodes(2)),'r--',x,Lessstd(:,Nodes(2)), 'r--', ...
    x,averages(:,Nodes(3)),'black-',x,Plusstd(:,Nodes(3)),'black--',x,Lessstd(:,Nodes(3)), 'black--')
legend(componentnames{1,Nodes(1)},'+/- std','', ...
    componentnames{1,Nodes(2)},'+/- std','', ...
    componentnames{1,Nodes(3)},'+/- std','')
xlim([0,x_max])
title({'Time evolution of average activity profile','across 100 repetitions for selected entities'})
xlabel('Timesteps')
ylabel('Activity profile')


figure(3)
x = 1:x_max;
Nodes=[10,20,34];
plot(x,averages(:,Nodes(1)),'b-',x,Plusstd(:,Nodes(1)),'b--',x,Lessstd(:,Nodes(1)), 'b--', ...
    x,averages(:,Nodes(2)),'r-',x,Plusstd(:,Nodes(2)),'r--',x,Lessstd(:,Nodes(2)), 'r--', ...
    x,averages(:,Nodes(3)),'black-',x,Plusstd(:,Nodes(3)),'black--',x,Lessstd(:,Nodes(3)), 'black--')
legend(componentnames{1,Nodes(1)},'+/- std','', ...
    componentnames{1,Nodes(2)},'+/- std','', ...
    componentnames{1,Nodes(3)},'+/- std','')
xlim([0,x_max])
title({'Time evolution of average activity profile','across 100 repetitions for selected entities'})
xlabel('Timesteps')
ylabel('Activity profile')

%%
% find which variables still evolve after steps 1000 (i.e. after
% perturbation released)

moving_var_average =    find((averages(1003,:) - averages(1000,:)) ~= 0);


%Sizes = cellfun('size',timeseries,1);

mycolors= ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30"];

figure(4)
%s= 66;
series = timeseries{s}; % also check 32, 55,66,86
x= 1:Sizes(s);

moving_var = find(sum(diff(series(1000:end,:,1)) ~= 0)~= 0);
moving_var(moving_var==12)=[];
half = ceil(length(moving_var)/2);

%plot( x([end-15:end]) , series( [end-15:end] , moving_var(1:half) , 1))
plot( x , series(:, moving_var(1:half) , 1))
legend(componentnames{1,moving_var(1:half)})
title({'Trajectories of variables still evolving after stopping', 'perturbation (i.e. after timestep #1000)'})
xlabel('Timesteps')
ylabel('Activity profile')
ylim([-0.2,1.2])%xline(1000)
xlim([995,x(end)])

% Set the line style order and color order
 ax= gca;
 ax.LineStyleCyclingMethod = "aftercolor";
 colororder(mycolors)
 ax.LineStyleOrder = ["-","--","-o"];


figure(5)
plot( x([end-20:end]) , series( [end-20:end] , moving_var(half+1:end) , 1));
legend(componentnames{1,moving_var(half+1:end)})
title({'Trajectories of variables still evolving after stopping ', 'perturbation (i.e. after timestep #1000)'})
xlabel('Timesteps')
ylabel('Activity profile')
ylim([-0.2,1.2])
%xline(1000)
ax= gca;
ax.LineStyleCyclingMethod = "aftercolor";
colororder(mycolors)
ax.LineStyleOrder = ["-","--","-o", "x-"];
xlim([995,x(end)])

%% Find which perturbation has SOX9 still evolving after perturb released

moving_SOX9 = [];
for s = 1:100
    serie = timeseries{s};
    if (serie(end,10,1)-serie(end-1,10,1)~=0)
        moving_SOX9 = [moving_SOX9, s];
    end

end

% result= SOX9 varies in perturbation 63 and 89 but 