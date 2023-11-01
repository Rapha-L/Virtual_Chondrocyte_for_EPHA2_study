function [attractors_test] = format_attractors(attractor_file)
% Utilities - funciton to format attractor profile in readable way for algorithm
    attractorsload = load(attractor_file,'attractors_read'); 
    attractors_test = attractorsload.attractors_read;
    %% if nodenames
    for i=1:length(attractors_test)
        attractors_test{i}=attractors_test{i}(2:4,:);
    end
    attractors_test = cellfun(@str2double,attractors_test,'UniformOutput',false);
    
    %save('attractorAC_WTdouble.mat','attractors')
end