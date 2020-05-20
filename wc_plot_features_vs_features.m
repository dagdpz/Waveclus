function handles=wc_plot_features_vs_features(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*5; %to prevent large plottings
if ~isfield('hfeatures',handles),
    if ~isempty(handles.features),
        handles=wc_create_featuresfig(handles);
        % Update handles structure
        if  ishandle(handles.mainfig)
            guidata(handles.mainfig, handles);
        end
    else
        return;
    end
end


interv=find(cellfun(@isempty,handles.classind) == 0);
featureind = [];
colind = [];
zeroind = max(interv);

for k=interv
    max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k}));
    featureind = cat(2,featureind,handles.classind{k}(randsample(length(handles.classind{k}),max_spikes))); % random sample... is actually important
    if k == zeroind && handles.ncl<k
        colind = cat(1,colind,zeros(max_spikes,3));
    else
        [X,~] = meshgrid(1:3,1:max_spikes);
        col=handles.colors(k,:);
        colind = cat(1,colind,col(X));
    end
end

colind = [colind ones(size(colind,1),1) * 0.5];

figure(handles.hfeatures);
used_colors=ismember(handles.colors,colind(:,1:3),'rows');
colormap(handles.colors(used_colors,:));
nf=length(handles.feature_names);
VER=version('-release');
VER=str2double(VER(1:4));
kk=2;
for i=2%:nf
    for j=i+1:nf
        cax=handles.hfaxes(kk);
        cla(cax);hold(cax,'on');
        %         hLine = plot(cax,handles.features(featureind,i),handles.features(featureind,j),'.','markersize',5);
        if VER>=2014 %since this part only works for matlab 2014 +
            hLine = scatter(cax,handles.features(featureind,i),handles.features(featureind,j),1,colind(:,1:3),'o');
            hLine.MarkerHandle.get;
            drawnow
            hLine.MarkerHandle.EdgeColorBinding = 'discrete';
            hLine.MarkerHandle.EdgeColorData = uint8(255*colind)';
        else
            [~, colix]=ismember(colind(:,1:3),handles.colors(used_colors,:),'rows');
            %colix=size(handles.colors,1)+1-colix;
            hLine = scatter(cax,handles.features(featureind,i),handles.features(featureind,j),30,colix,'.');
            drawnow
        end
        %         set(hLine.MarkerHandle,'EdgeColorBinding','discrete','EdgeColorData',uint8(255*colind)')
        clear hLine
        kk=kk+1;
    end
    
end