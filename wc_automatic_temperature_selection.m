function [handles,classind]=wc_automatic_temperature_selection(handles,classind,classtemp,tree)

min_clus=handles.WC.min_clus;
%% different algorithm:

% find peaks first

treetemp=cellfun(@numel,classtemp);
aa=diff(treetemp);
bb=diff((aa>0) + (aa<0)*-1)<0;
[ix1,ix2]=ind2sub(size(bb),find(bb));
ix1=ix1+1;


% detect valid peaks
clussiz=treetemp(sub2ind(size(classtemp),ix1,ix2));
val=clussiz>min_clus;
ix1=ix1(val);
ix2=ix2(val);
clussiz=clussiz(val);

% order peaks by size (smalles first) 

[~,sortix]=sort(clussiz);
ix1=ix1(sortix);
ix2=ix2(sortix);

% check how many bigger clusters without peaks there are and leave sapce for big clusters 
n_restclusters=min(ix2)-1;

if numel(ix1)>handles.WC.max_nrclasses-n_restclusters;    
ix1=ix1(numel(ix1)-handles.WC.max_nrclasses+n_restclusters+1:end);
ix2=ix2(numel(ix2)-handles.WC.max_nrclasses+n_restclusters+1:end);
end

% add restclusters
for c=1:n_restclusters
    ix1=[ix1; max(ix1)];
    ix2=[ix2; c];
end
% cluster all if no relevant cluster is found
if isempty(n_restclusters)&& isempty(ix1)
    ix1=1;
    ix2=1;
end


    n_classes=sum(~cellfun(@isempty,classind));
% take these as clusters and remove from the unsorted ones
for c=1:numel(ix1)
        t=classtemp{ix1(c),ix2(c)};
        t = setdiff(t,[classind{:}]); % only the ones not assigned yet
        n_classes=n_classes+1;
        classind{n_classes}=t;
    
end


classind(cellfun(@isempty,classind))=[];

%zero cluster, includes all unclustered spikes
nspk=numel(classtemp{1});
classind{end+1}=setdiff(1:nspk, [classind{:}]);

handles.WC.clus_per_temp=[ix1 ix2]';

% 
% n_classes=0;
% temp_start=1;
% handles.WC.clus_per_temp=[];
% while n_classes < handles.WC.max_nrclasses-1 && temp_start<size(tree,1)% max clusters ( leave 1 for unclustered) not reached yet
%     [temp] = wc_find_temperature(tree(temp_start:end,:),handles); %% need to reduce tree?
%     temp=temp+temp_start-1;
%     
%     %DEFINE CLUSTERS for specific temperature using min_clus variable
%     n_classes=sum(~cellfun(@isempty,classind));
%     clusters_for_this_temp=[];
%     for i=2:handles.WC.max_nrclasses,
%         t=classtemp{temp,i};
%         t = setdiff(t,[classind{:}]);
%         if length(t)>min_clus && n_classes < handles.WC.max_nrclasses
%             n_classes=n_classes+1;
%             classind{n_classes}=t;
%             clusters_for_this_temp=[clusters_for_this_temp i];
%         end
%     end
%     if temp==temp_start % didnt find appropriate temperature
%         break
%     end
%     if clusters_for_this_temp
%         handles.WC.clus_per_temp=[handles.WC.clus_per_temp [repmat(temp,1,numel(clusters_for_this_temp)); clusters_for_this_temp]];
%     end
%     temp_start=temp;
% end
% 
% 
% %% add biggest cluster of last iteration
% t=classtemp{temp,1};
% t = setdiff(t,[classind{:}]);
% if n_classes+1<handles.WC.max_nrclasses;
%     classind{n_classes+1}=t;
%     handles.WC.clus_per_temp=[handles.WC.clus_per_temp [temp; 1]];
% end
% 
% classind(cellfun(@isempty,classind))=[];
% 
% %zero cluster, includes all unclustered spikes
% classind{end+1}=setdiff(1:nspk, [classind{:}]);

end

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
end