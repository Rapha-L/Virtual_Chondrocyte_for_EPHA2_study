function []=Plot_bars_fun_attr(attractors,newdir)
% Create graph with computed wild type attractor profiles

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

string = ['component.mat'];
load(string) %same info also in screening_sup.attractors_read{i}(1,:)
xdata = {componentnames{1,:}};

basal = 1e-2; % for visualization with non nul bars

Fields= ["Healthy","Hypertrophic"];
labels= cellstr(Fields);

health =attractors{1}(1,:);
hyp = attractors{2}(1,:);



infla_nodes = [62,51,54,55,29];
ecm_nodes = [10,20,57,13,24,9];
    
Infl = [health(infla_nodes);hyp(infla_nodes)];
ecm = [health(ecm_nodes);hyp(ecm_nodes)];

Infl =Infl + basal;
ecm  = ecm + basal;

%ylim([0 1])

figure(1)
%bar(1:length(infla_nodes),Infl)
bar(Infl')
str = [xdata{infla_nodes}];
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ax = gca;
% ax.XTickLabel{1}=['\color{blue}' ax.XTickLabel{1}];
% ax.XTickLabel{2}=['\color{blue}' ax.XTickLabel{2}];
legend(labels)
fig1= gca;
saveas(fig1,fullfile(newdir,'Infl_nodes.fig'))

%ylim([0 1])

figure(2)
%bar(1:length(ecm_nodes),ecm)
bar(ecm')
str = [xdata{ecm_nodes}];
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ax = gca;
% ax.XTickLabel{1}=['\color{blue}' ax.XTickLabel{1}];
% ax.XTickLabel{2}=['\color{blue}' ax.XTickLabel{2}];
legend(labels)
fig3 = gca;
saveas(fig3,fullfile(newdir,'ECM_nodes.fig'))

end
