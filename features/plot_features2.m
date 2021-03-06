function handles=plot_features2(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*5; %to prevent large plottings
if ~isfield('hfeatures',handles),
    if ~isempty(handles.features),
        handles=create_featuresfig2(handles);
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
if length(interv) > 1
    interv = [interv(end) interv(1:end-1)];
end


for k=interv
    max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k}));
    featureind = cat(2,featureind,handles.classind{k}(1:max_spikes));
    if k == zeroind
        colind = cat(1,colind,zeros(max_spikes,3));
    else
        [X,~] = meshgrid(1:3,1:max_spikes);
        col=handles.colors(k,:);
        colind = cat(1,colind,col(X));
    end
end

colind = [colind ones(size(colind,1),1) * 0.1];

figure(handles.hfeatures);
colormap(handles.colors);
nf=length(handles.feature_names);
VER=version('-release');
VER=str2double(VER(1:4));
kk=1;
for i=1%:nf
    for j=i+1:nf
        cax=handles.hfaxes(kk);
        cla(cax);hold(cax,'on');
        %         hLine = plot(cax,handles.features(featureind,i),handles.features(featureind,j),'.','markersize',5);
        if VER>=2014 %since this part only works for matlab 2014 +
            hLine = scatter(cax,handles.features(featureind,i),handles.features(featureind,j),30,colind(:,1:3),'.');
            hLine.MarkerHandle.get;
            drawnow
            hLine.MarkerHandle.EdgeColorBinding = 'discrete';
            hLine.MarkerHandle.EdgeColorData = uint8(255*colind)';
        else
            [~, colix]=ismember(colind(:,1:3),flipud(handles.colors),'rows');
            colix=size(handles.colors,1)+1-colix;
            hLine = scatter(cax,handles.features(featureind,i),handles.features(featureind,j),30,colix,'.');
            drawnow
        end
        %         set(hLine.MarkerHandle,'EdgeColorBinding','discrete','EdgeColorData',uint8(255*colind)')
        clear hLine
        kk=kk+1;
    end
    
end