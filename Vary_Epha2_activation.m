%% Script for supplementary figure: use gradually increasing level of EPHA2 activation and study the effect on chondrocyte phenotype.

% This file is part of the repository Virtual_Chondrocyte_EPHA2_study
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2023 - KU Leuven

%File author(s): Raphaëlle Lesage

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
attractor_file = "attractorAC_WT_EPHA2.mat";

%% Define initial state
% From SOX9
attr_init= 1; %sox9+


attractorsload = load(attractor_file,'attractors'); 
All_attractors = attractorsload.attractors;
attr_init_profile = All_attractors{attr_init}; %pick initial state

filename = fullfile(newdir,'VaryEpha2_results.xls');   

%% Create structure with the in silico conditions to test
% Vary Epha2 activation with inflammation set to 1 (full activation)
condi = 1;
for Epha2= 0:0.1:1 %0:0.2:1
    for inflammation =1  % in case one also want to vary inflammation
        condition(condi).nodes = [62,51,54]; 
        condition(condi).inputvalues = [Epha2,inflammation,inflammation];
        condition(condi).initial_attr=attr_init ; %SOX9+
        condi = condi+1;
    end
end

name = repmat({'EPHA2_'},1,11);
inputvalue = string(1:11);  % string(0:0.1:1); but naming of field structure can't take '.' symbol replace by condition index (11 in total)
Fields= name+inputvalue; % brackets for necessary for plot funtion

for c = 1 : length(condition)
    field = Fields{c}
    %results_runs = Save_insilico_conditions_with_stats(condition(c).nodes,condition(c).inputvalues,condition(c).initial_attr,c,filename);
    [Transitionto_initAttrs,storeNewattr] = inSilicoConditions_return2(2/3,attractor_file,condition(c).initial_attr,condition(c).nodes,condition(c).inputvalues);
    condition(c).result_transition = Transitionto_initAttrs;
    condition(c).New_attractor = storeNewattr;
    [average_profile.(field),standDev.(field),w.(field)] = Average_output_profile(storeNewattr,Transitionto_initAttrs,All_attractors);

end

%% Compute average global activities from averge expression and average activation level

% this part could be integrated in the computation of average profile
% i.e.in Average_output_profile() (before creation of Average_destinationAttr)

average= average_profile;
for field = Fields
    average.(field)(1,:) = average.(field)(2,:).*average.(field)(3,:);
end 

%% Save output
file = fullfile(newdir,['VaryEpha2_experiment_results.mat']);
save(file,'condition','average_profile','standDev','w','average')

%Plot_bars_fun(attr_init_profile,average,standDev,newdir)

%Save_pie_chart(condition,Fields)


L= length(condition);
sox9 = zeros(1,L);
runx2 = zeros(1,L);

for c= 1: L
    sox9(c) = condition(c).result_transition(1);
    runx2(c) =  condition(c).result_transition(2);

end

epah2= 0:0.1:1;

plot(epah2,sox9,'-o', epah2, runx2, '-o')
ylim([0,110])
xlabel('EPHA2 activation level (input)')
ylabel({'Percentage of final states in each', 'chondrocyte profile'})
legend('Healthy state', 'Hypertrophic state')
title({'Effect of varying EPHA2 activation,','in combination with inflammation, on state transition'})

%% required functions
function[average_destinationAttr,standDev,warning]=Average_output_profile(storeNewattr,transition_percentage,attractors)
    %% average all the 'New attractors' to get the average value of each variable  
    warning = cell(1);
 if (transition_percentage==zeros(1,length(transition_percentage))) % all transitions end in new attractor
     for i= 1:length(storeNewattr)
        mat(1,:,i)= storeNewattr{i}(1,:);
        mat(2,:,i)= storeNewattr{i}(2,:);
        mat(3,:,i)= storeNewattr{i}(3,:);
     end
     
     average_globalFunction = mean(mat(1,:,:),3);
     average_proteinActivation = mean(mat(2,:,:),3);
     average_geneExpression = mean(mat(3,:,:),3);
     average_destinationAttr = [average_globalFunction;average_proteinActivation;average_geneExpression];
     standDev= std(mat,0,3);
     warning{1}= 'Only new attractors';
     
 elseif  (sum(transition_percentage) == 100)% all transitions go to one or several of the existing initial states
     a = transition_percentage ~=0;
     if size(find(a),2) >1  % more than one init_attractor as destination
        warning{1}='only init attractors as destination (more than one)';
         for i= 1:length(find(a))
             f = find(a);
            mat(1,:,i)= double(attractors{f(i)}(1,:));
            mat(2,:,i)= double(attractors{f(i)}(2,:));
            mat(3,:,i)= double(attractors{f(i)}(3,:));
         end

         W = transition_percentage(a);
         %average_globalFunction = mean(mat(1,:,:),3);
         average_globalFunction=compute_weighted_mean(mat(1,:,:),W);
         %average_proteinActivation = mean(mat(2,:,:),3);
         average_proteinActivation = compute_weighted_mean(mat(2,:,:),W);
         %average_geneExpression = mean(mat(3,:,:),3);
         average_geneExpression = compute_weighted_mean(mat(3,:,:),W);
         average_destinationAttr = [average_globalFunction;average_proteinActivation;average_geneExpression];
         standDev = std(mat,W,3);
         
     else % only one init_attractor as destination
         temp = double(attractors{a}(1:3,:));
         %W = transition_percentage;
         %temp = sum(temp.*W',1)/sum(W) 
         average_destinationAttr =temp;
         standDev = zeros(size(attractors{a}));  %take dimension of first attractor arbitrarily (all the same dimensions) 
         warning{1}= 'Only one of the initial attractors as destination';
     end
 else %mixture of init attractors and new attractors
     %cellfun(@isempty,storeNewattr)
      newAttr = find(~cellfun(@isempty,storeNewattr));
     for i= 1:length(newAttr)
        mat(1,:,i)= storeNewattr{newAttr(i)}(1,:);
        mat(2,:,i)= storeNewattr{newAttr(i)}(2,:);
        mat(3,:,i)= storeNewattr{newAttr(i)}(3,:);

     end
     
     start= length(newAttr);
     a = transition_percentage ~=0;
     f = find(a);
     for i= 1:length(f)       
            mat(1,:,start+i)= double(attractors{f(i)}(1,:));
            mat(2,:,start+i)= double(attractors{f(i)}(2,:));
            mat(3,:,start+i)= double(attractors{f(i)}(3,:));
     end
         W = [ones(1,length(newAttr)), transition_percentage(a)];
         average_globalFunction=compute_weighted_mean(mat(1,:,:),W);
         average_proteinActivation = compute_weighted_mean(mat(2,:,:),W);
         average_geneExpression = compute_weighted_mean(mat(3,:,:),W);
         average_destinationAttr= [average_globalFunction;average_proteinActivation;average_geneExpression];
         standDev = std(mat,W,3);
     warning{1}= 'Both New attractors and init attractors as destination';
 end


end
function [weighted_mean]=compute_weighted_mean(mat_row,W)
        for s =1:length(W)
            mat_row(1,:,s) = mat_row(1,:,s)*W(s);
        end
            Mysum = sum(mat_row(1,:,:),3);
        %(mat(1,:,1)*W(1) + mat(1,:,2)*W(2) + mat(1,:,3)*W(3))/sum(W);
        weighted_mean = Mysum/sum(W);
 end
