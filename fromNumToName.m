function[VarName]=fromNumToName(VarNum,ref)
% Varnum : int array
% VarName: string array
%ref: string providing the name of the matrice to be loaded (ex: 'componentnames.mat'
%or 'componentnames_AC.mat'

%This function translates variable numbers into their biochemical names

%%License
%This file is part of the repository Virtual_Chondrocyte_EPHA2_study
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

%load the list of component where idx is the varibale number
load(ref,'componentnames'); %structure with string array. eg. ref= 'compNames.mat'
componentnames=string(componentnames(1,:));
ntot=length(componentnames);
sz=size(VarNum);
%initialisation of the output string array
VarName=strings(sz);

%browse of the input varibale values
for i=1:sz(1)
    for j=1:sz(2)
        %make sure the input variable is within the scope of number of
        %existing component of the model (ntot=46 components)
        if (VarNum(i,j)>0 && VarNum(i,j)<=ntot) 
            %give the component name corresponding to input variable number
            VarName(i,j)=componentnames(VarNum(i,j)); 
        elseif (VarNum(i,j)>ntot+1 && VarNum(i,j)<=ntot*2) 
            %give the component name corresponding to input variable number
            VarName(i,j)=componentnames(VarNum(i,j)-46); 
        else
            VarName(i,j)=0; % set 0 if the input is not valid
        end
    end
end

end