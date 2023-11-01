%Script to analyze the result of the timeseries simulation after activation of EPHA2 and inflammation

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

%% Load results for EPHA2 activation plus inflammation
load('inSilico_condition (62  51)_1.mat')
load('component.mat')

% % find case with transition
% series=timeseries;
% for i=1:100
%     if series{i}(end,9,1) == 1
%         disp(i)
%     end
% end


Sizes = cellfun('size',timeseries,1);
max(Sizes)
min(Sizes)

x_max = min(Sizes);
storeMatrice = zeros(x_max,62,3,100);

for s =1:100
storeMatrice(:,:,:,s)=timeseries{s}(1:x_max,:,:);
end

%%
GlobalOnlyMatrice=zeros(max(Sizes),62,100);
GlobalOnlyMatrice(:,:,:)= storeMatrice(:,:,1,:);

averages = mean(GlobalOnlyMatrice,3);
stds = std(GlobalOnlyMatrice,0,3);
Plusstd = averages + stds;
Lessstd = averages - stds;


% find which variables still evolve after steps 1000 (i.e. after
% perturbation released)

moving_var =    find((averages(1003,:) - averages(1000,:)) ~= 0);

%% plot average of those variables
for f=1:length(moving_var)
figure(f)
x = 1:x_max;
plot( x(end-10:end), averages(end-10:end,f) )
ylim([0,1])
legend( componentnames{ 1 , moving_var(f) } )
end

%% get the full timeseries average for the moving variables


for m=moving_var
    s=100

    function [X,Y] = gettimeseries(timeseries,m,s)
        x_max = Sizes(s);
        Storematrice(:,:,:) =  timeseries{s}(1:x_max,m,:);
        GlobalOnlyMatrice(:,:)= storeMatrice(:,1);
    end
end
