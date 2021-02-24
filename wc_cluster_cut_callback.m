function wc_cluster_cut_callback(source, ~)
% added by Andrej Filippow
% open new window and allow to select lines to be removed from a cluster or
% points removed or added to the feature cluster

fig = figure(100);
cla;
copyobj(get(source,'children'),gca);

%% only to take over colors correctly
handles = guidata(get(source,'UserData'));
interv=find(cellfun(@isempty,handles.classind) == 0);
restcluster_present=handles.ncl<max(interv);
if restcluster_present
    cols=[handles.colors(interv(1:end-1),:) ; handles.colors(end,:)];
else
    cols=handles.colors(interv,:);
end
colormap(cols);
%set(0,'DefaultAxesColorOrder',handles.colors);


deletepoints = true;
h = struct;
Tag=get(source,'Tag');
set(fig,'Tag',Tag);
if strcmp(Tag,'Clus') % if called from cluster window
    set(fig,'Name','Select Lines to Delete');
end
if any(strfind(Tag,'Time_')) || any(strfind(Tag,'Feature1vs')) % if called from spike plot
    SelectionType=get(get(source,'Parent'),'SelectionType');
    if strcmp(SelectionType,'alt') % see if we want to create a new cluster ...
    set(fig,'Name','Delete Points');
        deletepoints = false;
    else % ...or remove from one
    set(fig,'Name','Select new Cluster');
    end
end

VER=version('-release');
VER=str2double(VER(1:4));
if VER>=2014
    h.featureIndex = [fig.Children.Children.UserData]; % tell which features it was plotting
else
    h.featureIndex = get(get(get(fig,'Children'),'Children'),'UserData');%
end
h.polygon = impoly();
if isvalid(h.polygon)
    set(fig,'CloseRequestFcn',{@closeCluster_callback h source deletepoints})
end

function closeCluster_callback(source, ~, h, outersource, deletepoints)
%return if user specified no points
if  ~isvalid(h.polygon)
    delete(source);
    return;
end
try
    %findobj(0, 'tag', 'GUI1');
    handles=guidata(get(outersource,'UserData'));
    nodes = getPosition(h.polygon);
    %return if the polygon is 2-Dimensional
    if length(nodes) < 3
        delete(source);
        return;
    end
%     % get rid of empty clusters, why were those allowed to exist in the
%     % first place?
%     validclusters = logical(cellfun(@numel,handles.classind));
%     validclusters(end) = 1;
%     handles.classind = handles.classind(validclusters);
    
    % next, get all clusters that have a point inside of the polygon
    if h.featureIndex{1} == 0 % cut spikes from spike plot
        ydata = handles.spikes(handles.classind{h.featureIndex{2}},:);
        xdata = repmat(-40+1:size(ydata,2)-40,size(ydata,1),1);
        indexdata = repmat(handles.classind{h.featureIndex{2}}', 1, size(ydata,2));
    else % add or remove from a feature plot
        xdata = [handles.features([handles.classind{:}], h.featureIndex{1})];
        if h.featureIndex{2}==1 % quick workaround for time
        ydata = [handles.index([handles.classind{:}])];
        else
        ydata = [handles.features([handles.classind{:}], h.featureIndex{2})];
        end
        indexdata = [handles.classind{:}];
    end
    
    index = inpolygon(xdata(:), ydata(:), nodes(:,1), nodes(:,2));
    hold on;
    scatter(xdata(index), ydata(index));
    drawnow;
    selectedPoints = unique(indexdata(index));
    
    % keep all clusters that are not empty after culling
    keepcluster = true(length(handles.classind),1);

    % rebuild clusters now
    if deletepoints
        for i = 1:length(handles.classind)-1
            handles.classind{i} = setdiff(handles.classind{i}, selectedPoints);
            if numel(handles.classind{i}) == 0
                keepcluster(i) = false;
                handles.ncl = handles.ncl - 1;
            end
        end
        handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);
        handles.classind = handles.classind(keepcluster);    
    else
        selectedPoints = setdiff(selectedPoints, [handles.classind{1:end-1}]);
        handles.classind = [selectedPoints handles.classind];
        handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);
        handles.ncl = handles.ncl + 1;
    end
    
    
 
    handles=wc_plot_spikes_and_ISI(handles);
    guidata(handles.mainfig, handles);
catch ex
    disp(ex.message)
end
set(handles.hclassify,'String','Classify');
set(handles.hclassify,'Value',0);
delete(source);
