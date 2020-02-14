function handles=plot_features(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*5; %to prevent large plottings
if ~isfield('hfeatures',handles),
    if ~isempty(handles.features), 
        handles=create_featuresfig(handles); 
        % Update handles structure
        guidata(handles.mainfig, handles);
    else
        return; 
    end
end
figure(handles.hfeatures);

nspk=handles.nspk;
ncl=handles.ncl;
nf=length(handles.feature_names);
kk=1;
for i=1:nf
    for j=i+1:nf
        cax=handles.hfaxes(kk);
        cla(cax);hold(cax,'on');
%         axes(handles.hfaxes(kk));cla;hold on
        if length(handles.classind{end})>nspk/2, interv=[0 1:ncl]; else interv=[1:ncl 0]; end
        for k=interv
            k1=mod(k-1,ncl+1)+1;
            max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k1}));
            if k==0, col=[0 0 0];
            else col=handles.colors(k,:);
            end
            x=handles.features(handles.classind{k1}(1:max_spikes),i);
            y=handles.features(handles.classind{k1}(1:max_spikes),j);
            if ~isempty(x)
%             plot(cax,x,y,['o' col],'MarkerFaceColor',[col 0.02],'markersize',1);
            hLine = plot(cax,x,y,'.','MarkerEdgeColor',col,'markersize',5);
            hLine.MarkerHandle.get;
            drawnow
            hLine.MarkerHandle.EdgeColorData = uint8(255*[col 0.1])';
            clear hLine 
            end
        end
        kk=kk+1;
        
    end
end
