function [d, sizes, da, percent, attractors] = Attractor_AC(nstates, saturation,input, inputindex)
% 1 argument + 3 optional arguments
% Function to compute the attractors of the systems (i.e. expression & activity profiles) and their canalization.

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

c = 1; % activity = c*PTM*transcription/translation/degradation
if nargin == 1
    input = [];
    inputindex = [];
    saturation = 2/3;
end
if nargin == 2
        input = [];
    inputindex = [];
end

if ~isnumeric(nstates)
   nstates = str2double(nstates);
end
%tic
% rng('shuffle')

global n
n = 62;
global tolerance
tolerance = 1e-3;


% saturation = 1;
global s
s = zeros(1,5);
s(1) = 1;
for i = 2:5
s(i) = 2* saturation/i;
end

slowidx = [5,11,13,16,20,21,23,24,25,44,46,47,48,49,51,59,61]; 
fastidx = [2,3, 7, 14, 18, 19, 22, 26, 29, 30, 35, 36, 40, 41, 45, 50, 54, 55, 56];

d = zeros(1,nstates); 

k = 1;  % k is the number of the attractor
ptr = 1; % counts the total amount of initial states we have assessed

attractors = {ones(35,1)}; % will be used to store the attractors
attrcnt = 0;
aSox9 = 0;
aRunx2 = 0;
aMMP13= 0;
%wb = waitbar(0,'computing attractors and distances...');
while ptr <= nstates
    ptr
    %waitbar(ptr/nstates,wb);

    newstate = zeros(3,n); 
    newstate(1,:) = rand(n,1);
    newstate(1,inputindex) = input; % keep inputnodes fixed
    newstate(2,:) = sqrt(newstate(1,:));
    newstate(3,:) = newstate(2,:);
    newstate(2,slowidx) = 1; % full control of fast controlled nodes in fast variable, same for slow
    newstate(3,slowidx) = newstate(1,slowidx);
    newstate(3,fastidx) = 1;
    newstate(2,fastidx) = newstate(1,fastidx);

    attr = 0;
    cnt = 1; % indicates the number of the state in the new sequence (reinitialises after every new random state)
 while (attr == 0)
    attr = 1;
    for l = randperm(n)
        attrf = 0;
%     hh = 0;
    while (attrf == 0) % && hh < 100
        % hh = hh + 1;
        attrf = 1;
        for j = randperm(n)
            changednode = zupdateas_AC_EPHA2(newstate,j,'fast'); 
        if (abs(newstate(2,j)-changednode) > tolerance) && ~ismember(j,inputindex) % newstate(2,j) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
        
           newstate(2,j) = changednode;
           newstate(1,j) =  min(1,c*newstate(2,j)*newstate(3,j)); % min(newstate(2,j), newstate(3,j));
           attrf = 0;
           newstate(1,inputindex) = input; % keep inputnodes fixed
           
           break;
        end
        end
    end
            changednode = zupdateas_AC_EPHA2(newstate,l,'slow');
            %if new state different from previous
         if (abs(newstate(3,l)-changednode) > tolerance) && ~ismember(l,inputindex) %newstate(3,l) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
           newstate(3,l) = changednode;
           newstate(1,l) = min(1,c*newstate(2,l)*newstate(3,l));
           attr = 0;
           cnt = cnt + 1;
           newstate(1,inputindex) = input; % keep inputnodes fixed
           
           break;
         end
    end     
 end
   % we need to 'remove' the last state, this was an attractor, by lowering
   % pointer we will overwrite it

   % 2 options: we encountered attractor before, in which case we need to
   % know which k it corresponds to
   if checkstates(newstate,attractors,k-1)
      [~,old] = checkstates(newstate,attractors,k); % we encountered attractor before
            d(ptr) = cnt;   % reverse the distances to the attractor
%             ab(ptr-cnt+1:ptr) = old;      % we go back and change everything upstream to match the old attractor
            attrcnt(old) = attrcnt(old) + 1;
      % -attrk(old), if every attractor is a singleton we can just use
      % index like above
   else % we found a new attractor
       d(ptr) = cnt;   % reverse the distances to the attractor
%        ab(ptr-cnt+1:ptr) = k;      % we go back and change everything upstream to match the old attractor
       attractors(k) = {newstate};
       attrcnt(k) = 1;
       aRunx2(k)=newstate(1,9);
       aSox9(k)=newstate(1,10);
       aMMP13(k)=newstate(1,24);
       k = k + 1;
   end
   ptr = ptr + 1;
end


%close(wb)

%% Determine size basins and values of Sox9/Runx2

% [~, attrcount] = count_unique(ab); 
k = k-1; % since we never found an attractor to match the last k
sizes = zeros(6,k);
sizes(1,1:k) = 1:k; % numbering the attractors
% sizes(2,:) = attrcount(k:-1:1); % the number of states in the attractor (ie singleton vs cyclic)
% sizes(3,:) = attrcount(1:k); % basin size
sizes(3,:) = attrcnt;
sizes(4,:) = aSox9; % Determine the amount of states in the attractor that are Sox9 positive: max is 1 * size attractor (1 sox9 in every state)
sizes(5,:) = aRunx2; 
sizes(6,:) = aMMP13;

% percentages in Runx2, Sox9, both, none
percent = cell(4,2);
percent(1,1) = {'Runx2'};
percent(2,1) = {'Sox9'};
percent(3,1) = {'both'};
percent(4,1) = {'none'};
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
for i = 1:k


if sizes(4,i) > 0.05
    if sizes(5,i) > 0.05
    temp3 = temp3 + sizes(3,i)/sum(sizes(3,1:k));
    else temp2 = temp2 + sizes(3,i)/sum(sizes(3,1:k));
    end
    else if sizes(5,i) > 0.05
            temp1 = temp1 + sizes(3,i)/sum(sizes(3,1:k));
        else temp4 = temp4 + sizes(3,i)/sum(sizes(3,1:k));
        end
end
end
    percent(1,2) = {temp1};
    percent(2,2) = {temp2};
    percent(3,2) = {temp3};
    percent(4,2) = {temp4};       
% average distance to attractor: algorithm won't count states that start on
% an attractor, but is only insignificant fraction of state space
da = zeros(size(d));
    da(1) = d(1);
for i= 1:length(d) - 1
    if d(i) <= d(i+1) % greater than or equal(if the new state is 1 removed from attractor it would'nt be counted otherwise)
        da(i) = d(i+1); % You will reach this when on an attractor (0 = 0) but will be removed later on anyway
    end
end
da(da < 1) = []; % remove all empty elements from da
da = mean(da);
string = ['attractorAC_WT_EPHA2.mat']; 
attractors_read=attractors;
% formatting for easy reading:
load('component.mat','componentnames')
for a = 1:length(attractors_read)
    attractors_read{a}=[componentnames{1,:}; attractors_read{a}];
end
save(string, 'd', 'sizes', 'da', 'percent', 'attractors','attractors_read');

%t = toc

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