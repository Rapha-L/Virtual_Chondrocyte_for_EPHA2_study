%% Analyze screening of EPHA2 in combi with a second entity

load('inSilico_EphA2_screen_initSOX9.mat')
% Analyze screening of perturbations from the Sox9+ attractor
temp_cell = struct2cell(condition_Sox9);

transition_results=cell(124,1);
transition_results(:,1) = temp_cell(5,1,:);


hit_1 =zeros(124,1);
hit_2 =zeros(124,1);

for i =1:length(transition_results)    
    if transition_results{i}(1)<80
        hit_1(i)=1;
        if transition_results{i}(1)>=2
            hit_2(i)=1; % remove 100% transtiion to make the graph lighter 
        end
    end
    
end

hit_1 = logical(hit_1);
PotentiateEPHA2 = condition_Sox9(hit_1);

% create lables
load('component.mat')
Labels =cell(124,1);
for L=1:length(componentnames)
    L_temp = upper(componentnames{1,L});
    Labels{L} = [char(L_temp) 'up'];
    Labels{L+length(componentnames)} = [char(L_temp) 'down'];
end

toto = 'Labels';
saved_Labels= Labels(hit_1);
[PotentiateEPHA2.(toto)]=saved_Labels{1:end};

Labels_filtered_1 = Labels;
Labels_filtered_2 = Labels;

index=find(hit_1);
index2=find(hit_2);

for i=1:length(Labels_filtered_1 )
    if ~ismember(i,index)
        Labels_filtered_1{i}=" ";
    end
    
    if ~ismember(i,index2)
        Labels_filtered_2{i}=" ";
    end
    
end

%% plot 

%S= struct
X = zeros(1,124);
Y = zeros(1,124);

for i =1:length(transition_results)    
   X(i) = transition_results{i}(1);
   Y(i) = transition_results{i}(2);
end
c = autumn; %linspace(1,10,length(X));
s=scatter(X,Y,10,'filled');
ylabel('Propability of transition towards Runx2+');
xlabel('Propability of remaining in Sox9+');
title({'Likehood for transition from the Sox9+ state' ,'under EPHA2 inhibition combined with a second perturbation'},'Fontsize',12) 

t = cellstr(Labels_filtered_2);
dx = 1; dy = 1; % displacement so the text does not overlay the data points
text(X+dx, Y+dy,t,'Fontsize', 9);

