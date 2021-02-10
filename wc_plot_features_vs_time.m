function handles=wc_plot_features_vs_time(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*10;
if ~isfield('htimecourse',handles) || ~ishandle(handles.htimecourse),
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
toplotind = [];
colind = [];
zeroind = max(interv);

total_N_spikes= numel([handles.classind{:}]);
for k=interv
%     max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k}));
    max_spikes=round(length(handles.classind{k})/total_N_spikes*MAX_SPIKES_TO_PLOT);
    spk_indexes=randperm(length(handles.classind{k}));
    toplotind = cat(2,toplotind,handles.classind{k}(spk_indexes(1:max_spikes)));    
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

%%first plot is histogram
cax=handles.htaxes(1);
cla(cax);hold(cax,'on');
mxx=prctile(handles.features(:,1),99.9)*1.1;
mnx=prctile(handles.features(:,1),0.1)*1.1;
bins=[min(handles.index(toplotind)) max(handles.index(toplotind))];
bins=[bins(1):diff(bins)/100:bins(2)-diff(bins)/200]+diff(bins)/200;

    histo=hist(handles.index,bins);
    max_all=max(histo);
    histo=histo/max_all*(mxx-mnx);
    histo=histo+mnx;
    plot(cax,histo,bins,'k:','linewidth',2);
for k=interv
    histo=hist(handles.index(handles.classind{k}),bins);
    histo=histo/max_all*(mxx-mnx);
    histo=histo+mnx;
    plot(cax,histo,bins,'color',handles.colors(mod(k-1,size(handles.colors,1))+1,:));
end


for i=2:nf
    %for j=i+1:nf
    cax=handles.htaxes(i);
    cla(cax);hold(cax,'on');
        
    %         hLine = plot(cax,handles.features(featureind,i),handles.features(featureind,j),'.','markersize',5);
    if VER>=2014 %since this part only works for matlab 2014 +
        hLine = scatter(cax,handles.features(toplotind,i),handles.index(toplotind),1,colind(:,1:3),'o','filled');
        hLine.MarkerHandle.get;
        drawnow
        hLine.MarkerHandle.EdgeColorBinding = 'discrete';
        hLine.MarkerHandle.EdgeColorData = uint8(255*colind)';
        hLine.MarkerHandle.FaceColorBinding = 'discrete';
        hLine.MarkerHandle.FaceColorData = uint8(255*colind)';
        hLine.UserData = {1 i};
    else
            [~, colix]=ismember(colind(:,1:3),handles.colors(used_colors,:),'rows');
        hLine = scatter(cax,handles.features(toplotind,i),handles.index(toplotind),30,colix,'.');
        set(hLine,'UserData',{1 i});
        drawnow
    end
    %         set(hLine.MarkerHandle,'EdgeColorBinding','discrete','EdgeColorData',uint8(255*colind)')
    clear hLine
    %kk=kk+1;
    %end
    
end