function wc_cluster_cut_callback(source, ~)
% added by Andrej Filippow
% open new window and allow to select lines to be removed from a cluster or
% points removed or added to the feature cluster

fig = figure(100);
cla;
copyobj(get(source,'children'),gca);

deletepoints = true;
h = struct;
SelectionType = 'none';
handles = guidata(get(source,'UserData'));

%% find out from where we calle dit
Tag=get(source,'Tag');
set(fig,'Tag',Tag);
if strfind(Tag,'Clus') % if called from cluster window
    set(fig,'Name','Select Lines to Delete');
    interv=[];    
    aa=get(get(fig,'Children'),'Children');    
    bb=get(aa,'UserData');
    h.featureIndex =bb{~cellfun(@isempty,bb)};
else  % if called from features figure
    %% only to take over colors correctly
    interv=find(cellfun(@isempty,handles.classind) == 0);
    restcluster_present=handles.ncl<max(interv);
    if restcluster_present
        cols=[handles.colors(interv(1:end-1),:) ; handles.colors(end,:)];
    else
        cols=handles.colors(interv,:);
    end
    colormap(cols);    
    stepsize=1/numel(interv);
    xpos=0;
    for c=interv
        handles.(['clus' num2str(c)])=uicontrol('units','normalized','Style','togglebutton','String',['clus' num2str(c)],'FontSize',12,...
            'Tag',num2str(c),'BackgroundColor',cols(c,:),...
            'Position',[stepsize*1/6+xpos 0.97 stepsize*2/3 0.035],...
            'Callback',{@clusterselection_Callback});
        xpos=xpos+stepsize;
    end    
    aa=get(get(fig,'Children'),'Children');
    h.featureIndex = get(aa{~cellfun(@isempty,aa)},'UserData');%
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

if strcmp(SelectionType,'extend') % see if we want to simply cut straight
    h.polygon = imrect();
else
    h.polygon = impoly();
end
if isvalid(h.polygon)
    set(fig,'CloseRequestFcn',{@closeCluster_callback h deletepoints handles interv})
end

function closeCluster_callback(source, ~, h, deletepoints, handles, interv)
%return if user specified no points
valid_clusters=zeros(size(handles.classind));
for c=interv
    valid_clusters(c)=get(handles.(['clus' num2str(c)]),'Value');
end
if ~any(valid_clusters)
    valid_clusters=ones(size(valid_clusters));
end
if  ~isvalid(h.polygon)
    delete(source);
    return;
end
try
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
        xdata = repmat(1:size(ydata,2),size(ydata,1),1);
        indexdata = repmat(handles.classind{h.featureIndex{2}}', 1, size(ydata,2));
        index = inpolygon(xdata(:), ydata(:), nodes(:,1), nodes(:,2));
    else % add or remove from a feature plot
        xdata = [handles.features([handles.classind{:}], h.featureIndex{1})];
        if h.featureIndex{2}==1 % quick workaround for time
            ydata = [handles.index([handles.classind{:}])];
        else
            ydata = [handles.features([handles.classind{:}], h.featureIndex{2})];
        end
        indexdata = [handles.classind{:}];
        index = inpolygon(xdata(:), ydata(:), nodes(:,1), nodes(:,2));
        [validindex]=ismember([handles.classind{:}],[handles.classind{valid_clusters==1}]);
        index = index & validindex';
    end
    
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

function clusterselection_Callback(source,~)
