function [temp] = wc_find_temperature(tree,handles)
% Selects the temperature.
min_clus=handles.WC.min_clus;
diff_matrix=diff(tree(:,5:end));
num_temp=size(diff_matrix,1);
% ORIGINAL
temp = 1;         % Initial value
for t=1:num_temp-1;
    if any(diff_matrix(t,:)>min_clus & diff_matrix(t+1,:)<min_clus) %% what if its the last one?
        temp=t+1;     
        break;
    end
end
   
