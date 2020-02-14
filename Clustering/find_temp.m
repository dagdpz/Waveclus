function [temp] = find_temp(tree,handles)
% Selects the temperature.

num_temp=handles.par.num_temp;
min_clus=handles.par.min_clus;

aux =diff(tree(:,5));   % Changes in the first cluster size
aux1=diff(tree(:,6));   % Changes in the second cluster size
aux2=diff(tree(:,7));   % Changes in the third cluster size
aux3=diff(tree(:,8));   % Changes in the third cluster size
temp = 1;               % Default value in case no more than one cluster appears.

for t=1:num_temp-1;
    % Looks for the appearance of a cluster larger than min_clus.
    if ( aux(t) > min_clus || aux1(t) > min_clus || aux2(t) > min_clus || aux3(t) >min_clus )    
        temp=t+1;         
    end
end
   

% % for t=1:num_temp-1;
% % % A threshold lager than min_clus in the first 2 clusters means that other cluster appeared.
% %     if ( aux(t) < -min_clus & (aux1(t) > min_clus | aux2(t) > min_clus) ) |...
% %        ( aux1(t) < -min_clus & (aux2(t) > min_clus))    
% % %        if tree(t+2,6) >= stab * tree(t+1,6)   %Also requires some stability of the 2nd cluster.
% %             temp=t+1;         
% % %            break           % Assumes that clusters appear at the same temperature.
% % %        end
% %     end
% % end
