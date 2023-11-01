function [perresult, stabattr, trlik,storeNewattr,transition,twicepresult,distribution] = StateStab(saturation,runnumber,attractor_file)
% STATESTAB Calculates the stability of a state by perturbing nodes one by
% one

%% License
%This file is part of the repository Virtual_Chondrocyte_EPHA2_study
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2023 - KU Leuven

%Distributed under the GPLv3 License.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>


%%

[attractors] = format_attractors(attractor_file);
disp(['formatting done...'])
runs = 100;
fixedtime = 1000;

if ~isnumeric(saturation)
   saturation = str2num(saturation);
end

global n
n = 62;   
global s
s = zeros(1,5);
s(1) = 1;
for i = 2:5
s(i) = 2* saturation/i;
end

perresult = cell(length(attractors),n,2); 
storeNewattr= cell(length(attractors),n,2);
twicepresult = cell(length(attractors),n,2); %store the result without dividing per 2 {attrcnt}for some nodes


tic
for i = 1:length(attractors)
    disp(['attractor ' num2str(i) ' started...'])
    attr = attractors{i}; 
% attractor{i} is a matrix of size=(3*46) where the first line is the node activity and each colone corresponds to one node
    for j = [1:n] 
        
        disp(['node ' num2str(j) ' ...'])
        if (attr(1,j) ~= 1) && (attr(1,j) ~= 0) % 2 possible flip states: 0 and 1: first 1 % rest algoritme uitwerken
        for k = 1:2
            per = attr;
            per(1,j) = k-1;% create perturbed state
            if ~isequal(attr,zupdates_AC_EPHA2(per))
                [attrcnt,~,stateruns] = ssattr(per, runs, attractors, fixedtime,j);
                perresult(i,j,k) = {attrcnt/2}; % store results of perturbation
                twicepresult(i,j,k) = {attrcnt};
                storeNewattr(i,j,k)={stateruns};  
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                twicepresult(i,j,k) = {attrcnt};
                attrcnt(i) = runs/2; % we want to give this the same weight: since every node has weight one, and we perform 2 simulations here, we divide by 2
                perresult(i,i,k) = {attrcnt};
                storeNewattr(i,j,k)={'no transition'};  
            end
            
        end
        else % we only explore 2 values if the initial one is intermediate
        per = attr;
        per(1,j) = ~per(1,j);% create perturbed state
            if ~isequal(attr,zupdates_AC_EPHA2(per))
                [attrcnt,~,stateruns] = ssattr(per, runs, attractors, fixedtime,j);
                twicepresult(i,j,int8(-attr(1,j)+2)) = {attrcnt};
                perresult(i,j,int8(-attr(1,j)+2)) = {attrcnt}; % store results of perturbation
                storeNewattr(i,j,int8(-attr(1,j)+2))={stateruns};
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                attrcnt(i) = runs; % we want to give this the same weight
                twicepresult(i,j,int8(-attr(1,j)+2)) = {attrcnt};
                perresult(i,i,int8(-attr(1,j)+2)) = {attrcnt};
                storeNewattr(i,j,int8(-attr(1,j)+2))={stateruns};
            end
        
        end
    
    end
    
end
% 


% stability per attractor: count all flipped nodes with effect and percentages
% perresultav = perresult(:,:,1) + perresult(:,:,2);% this is what we get if you flip intermediate values randomly
nattr = length(attractors); % determine the number of attractors
stabattr = cell(3,nattr,2); % stores how many flips have the potential to change attractors, second row contains overall prob, third row probability per state
trlik = zeros(nattr,nattr,2); % contains likelihood of transition (going from attr to another attr)
for c = 1:2
for a = 1:nattr
     flip = 0; % contains the number of flipped genes
     perflip = zeros(1, n); % contains the likelihood of flipping per gene
     for b = 1:n
        if ~isempty(twicepresult{a,b,c})
        x = twicepresult{a,b,c};
        if (x(a) ~= runs)
            flip = flip + 1;
            perflip(b) = (sum(x) - x(a))/runs; % percentage of states flipping
            trlik(a,1:nattr,c) = trlik(a,1:nattr,c) + x(1:nattr); % add number of transitions to transition matrix
        else
            trlik(a,a,c)= trlik(a,a,c) + runs; % all runs end up in the same and original attractor
        end
        end
     end
%      trlik(a,:) = trlik(a,:)/sum(trlik(a,:));
     stabattr(1,a,c) = {flip};
     stabattr(3,a,c) = {perflip};
     stabattr(2,a,c) = {sum(perflip)/(flip*runs)}; % average chance of flipping per node that can potentially switch states
end
end
toc
trlik2 = mean(trlik,3)*2;


%calculation of matrix of transition.
transition=trlik2./sum(trlik2,2); % sum(trlik2,2)not always equals n*runs(=4600)
distribution=stationary(transition);
%save variables in file.mat
string = ['single_perturbation' num2str(runnumber) 'sat' num2str(saturation) '.mat'];
save(string, 'stabattr', 'perresult', 'trlik', 'trlik2','saturation','storeNewattr','transition','distribution')
   

end

function [attrcnt,newattr,stateruns] = ssattr(initstate,runs,attractors, fixedtime, fixednode)
% SSATTR single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state

stateruns=cell(runs,1);
attrcnt = zeros(1,length(attractors));
global tolerance
global n
tolerance = 1e-2; % should be somewhat higher than attr finding algorithm
newattr=0;
for z = 1:runs

t = 0;
state = initstate;
attr = 0;  
 while (attr == 0)
    t = t + 1;

    if t > fixedtime
    attr = 1; % you can only reach attractor after perturbation period
    state(1,:) = state(2,:) .* state(3,:);
    end
    temp = randperm(n);
    if t <= fixedtime % depending on whether fixed interval passed or not
    temp(temp == fixednode) = []; % remove the node that should stay fixed
    end
        
    attrf = 0;
% if change = 1 % do fast loop only if there was a change in the slow one
    while (attrf == 0) 
        attrf = 1;
        tmp = randperm(n);
        if t <= fixedtime % depending on whether fixed interval passed or not
        tmp(tmp == fixednode) = []; % remove the node that should stay fixed
        end
        for j = tmp
            changednode = zupdateas_AC_EPHA2(state,j,'fast'); 
        if (abs(state(2,j)-changednode) > tolerance)  %&& ~ismember(j,inputindex) % newstate(2,j) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
           state(2,j) = changednode;
           state(1,j) =  state(2,j)*state(3,j); % min(newstate(2,j), newstate(3,j));
           attrf = 0;
%           state(1,inputindex) = input; % keep inputnodes fixed
           break;
        end
        end
    end
%  end
        
    for l = temp
            changednode = zupdateas_AC_EPHA2(state,l,'slow');
         if (abs(state(3,l)-changednode) > tolerance) %&& ~ismember(l,inputindex) %newstate(3,l) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
           state(3,l) = changednode;
           state(1,l) = state(2,l)*state(3,l);
           attr = 0;
%            change = 1;
           
            %          state(1,inputindex) = input; % keep inputnodes fixed
           break;
         end
         
    end

 end
     if ~checkstates(state,attractors,length(attractors))
        disp('new attractor found')
        newattr=newattr+1;
        stateruns{z}=state;
    else
        [~,old] = checkstates(state,attractors,length(attractors));
        attrcnt(old) = attrcnt(old) + 1;
    end  
end
end

function [boolean, number] = checkstates(newstate,state,maxidx)
% function to check whether a state has already been described
global tolerance
boolean = 0;
number = 0;

if maxidx ~= 0
    for i = 1:maxidx
        if  abs(newstate-state{i}) < 100*tolerance % isequal(newstate, state{i});
            boolean = 1;
            number = i;
            
            break;
        end
    end
end
end
