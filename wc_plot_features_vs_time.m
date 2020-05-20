function handles=wc_plot_features_vs_time(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*5; %to prevent large plottings
if ~isfield('htimecourse',handles),
    if ~isempty(handles.index),
        handles=wc_create_featurevstimefig(handles);
        % Update handles structure
        if  ishandle(handles.mainfig)
            guidata(handles.mainfig, handles);
        end
    else
        return;
    end
end


interv=find(cellfun(@isempty,handles.classind) == 0); %different clusters
featureind = [];
colind = [];
zeroind = max(interv);


for k=interv
    max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k}));
    spk_indexes=randperm(length(handles.classind{k}));
    featureind = cat(2,featureind,handles.classind{k}(spk_indexes(1:max_spikes)));
    if k == zeroind &&  handles.ncl<k
        colind = cat(1,colind,zeros(max_spikes,3));
    else
        [X,~] = meshgrid(1:3,1:max_spikes);
        col=handles.colors(mod(k-1,size(handles.colors,1))+1,:);
        colind = cat(1,colind,col(X));
    end
end

colind = [colind ones(size(colind,1),1) * 0.5];

figure(handles.htimecourse);
used_colors=ismember(handles.colors,colind(:,1:3),'rows');
colormap(handles.colors(used_colors,:));
nf=length(handles.feature_names);
VER=version('-release');
VER=str2double(VER(1:4));

%kk=1;
for i=1:nf
    %for j=i+1:nf
    cax=handles.htaxes(i);
    cla(cax);hold(cax,'on');
    %         hLine = plot(cax,handles.features(featureind,i),handles.features(featureind,j),'.','markersize',5);
    if VER>=2014 %since this part only works for matlab 2014 +
        hLine = scatter(cax,handles.features(featureind,i),handles.index(featureind),1,colind(:,1:3),'o','filled');
        hLine.MarkerHandle.get;
        drawnow
        hLine.MarkerHandle.EdgeColorBinding = 'discrete';
        hLine.MarkerHandle.EdgeColorData = uint8(255*colind)';
        hLine.MarkerHandle.FaceColorBinding = 'discrete';
        hLine.MarkerHandle.FaceColorData = uint8(255*colind)';
    else
            [~, colix]=ismember(colind(:,1:3),handles.colors(used_colors,:),'rows');
        hLine = scatter(cax,handles.features(featureind,i),handles.index(featureind),30,colix,'.');
        drawnow
    end
    %         set(hLine.MarkerHandle,'EdgeColorBinding','discrete','EdgeColorData',uint8(255*colind)')
    clear hLine
    %kk=kk+1;
    %end
    
end