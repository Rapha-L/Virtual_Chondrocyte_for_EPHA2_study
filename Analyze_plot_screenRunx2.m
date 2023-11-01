%% Analyze screening of EPHA2 in combi with perturbation of a second entity (from Runx2+)

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


%% Load results of screening from the hypertrophic chondrocyte
load('inSilico_EphA2-_screen_initRUNX2.mat')
n =122;
% Analyze screening of perturbations from the Runx2+ attractor
temp_cell = struct2cell(condition);

transition_results=cell(n,1);
transition_results(:,1) = temp_cell(5,1,:);


hit_1 =zeros(n,1);
% hit_2 =zeros(n,1);

for i =1:length(transition_results)    
    if transition_results{i}(1)>0
        hit_1(i)=1;
%         if transition_results{i}(1)>=2
%             hit_2(i)=1; % remove 100% transtiion to make the graph lighter 
%         end
    end
    
end

hit_1 = logical(hit_1);
PotentiateEPHA2 = condition(hit_1);

% create lables
load('component.mat')
Labels =cell(n,1);
for L=1:length(componentnames)-1
    L_temp = upper(componentnames{1,L});
    Labels{L} = [char(L_temp) 'up'];
    Labels{L+length(componentnames)-1} = [char(L_temp) 'down'];
end

toto = 'Labels';
saved_Labels= Labels(hit_1);
[PotentiateEPHA2.(toto)]=saved_Labels{1:end};

Labels_filtered_1 = Labels;
%Labels_filtered_2 = Labels;

index=find(hit_1);
%index2=find(hit_2);

for i=1:length(Labels_filtered_1 )
    if ~ismember(i,index)
        Labels_filtered_1{i}=" ";
    end
    
%     if ~ismember(i,index2)
%         Labels_filtered_2{i}=" ";
%     end
    
end

%% plot 

%S= struct
X = zeros(1,n);
Y = zeros(1,n);

for i =1:length(transition_results)    
   X(i) = transition_results{i}(2);
   Y(i) = transition_results{i}(1);
end
c = autumn; %linspace(1,10,length(X));
s=scatter(X,Y,10,'filled');
ylabel('Propability of transition towards Sox9+');
xlabel('Propability of remaining in Runx2+');
title({'Likehood for transition from the Runx2+ state' ,'under EPHA2 activation combined with a second perturbation'},'Fontsize',12) 

t = cellstr(Labels_filtered_1);
dx = 0.2; dy = -0.2; % displacement so the text does not overlay the data points
text(X+dx, Y+dy,t,'Fontsize', 9);

