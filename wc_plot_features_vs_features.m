function handles=wc_plot_features_vs_features(handles)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT*5; %to prevent large plottings
if ~isfield(handles,'hfeatures') %|| ~isvalid(handles.hfeatures),
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
toplotind = [];
colind = [];
zeroind = max(interv);

for k=interv
    max_spikes = min(MAX_SPIKES_TO_PLOT,length(handles.classind{k}));
    toplotind = cat(2,toplotind,handles.classind{k}(randsample(length(handles.classind{k}),max_spikes))); % random sample... is actually important
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

for i=2:3%:nf
    cax=handles.hfaxes(i-1);
    cla(cax);hold(cax,'on');
    y_lim=get(cax,'ylim');
    mny=y_lim(1);
    mxy=y_lim(2);
    x_lim=get(cax,'xlim');
    mnx=x_lim(1);
    mxx=x_lim(2);
    % mxx=prctile(handles.features(:,2),99.9)*1.1;
    % mnx=prctile(handles.features(:,2),0.1)*1.1;
    bins=[prctile(handles.features(toplotind,i),0.1)*1.1 prctile(handles.features(toplotind,i),99.9)*1.1];
    if i==2
        binfactor=(mxx-mnx)/diff(bins);
        binoffset=mnx-bins(1)*binfactor;
        mx=mxy;
        mn=mny;
    else
        mx=mxx;
        mn=mnx;
        binfactor=(mxy-mny)/diff(bins);
        binoffset=mny-bins(1)*binfactor;
        
    end
    bins=[bins(1):diff(bins)/100:bins(2)-diff(bins)/200]+diff(bins)/200;
    
    histo=hist([handles.features(:,i)],bins);
    max_all=max(histo);
    histo=histo/max_all*(mx-mn);
    histo=histo+mn;
        if i==2
    plot(cax,bins*binfactor+binoffset,histo,'k:','linewidth',2);
        else
    plot(cax,histo,bins*binfactor+binoffset,'k:','linewidth',2);
        end
    for k=interv
        histo=hist(handles.features(handles.classind{k},i),bins);
        histo=histo/max_all*(mx-mn);
        histo=histo+mn;
        if i==2
            plot(cax,bins*binfactor+binoffset,histo,'color',handles.colors(mod(k-1,size(handles.colors,1))+1,:));
        else
            plot(cax,histo,bins*binfactor+binoffset,'color',handles.colors(mod(k-1,size(handles.colors,1))+1,:));
        end
    end
end
kk=3;
for i=2%:nf
    for j=i+1:nf
        cax=handles.hfaxes(kk);
        cla(cax);hold(cax,'on');
        %         hLine = plot(cax,handles.features(featureind,i),handles.features(featureind,j),'.','markersize',5);
        if VER>=2014 %since this part only works for matlab 2014 +
            hLine = scatter(cax,handles.features(toplotind,i),handles.features(toplotind,j),1,colind(:,1:3),'o');
            hLine.MarkerHandle.get;
            drawnow
            hLine.MarkerHandle.EdgeColorBinding = 'discrete';
            hLine.MarkerHandle.EdgeColorData = uint8(255*colind)';
            hLine.UserData = {i j};
        else
            [~, colix]=ismember(colind(:,1:3),handles.colors(used_colors,:),'rows');
            %colix=size(handles.colors,1)+1-colix;
            hLine = scatter(cax,handles.features(toplotind,i),handles.features(toplotind,j),30,colix,'.');
            set(hLine,'UserData',{i j});
            drawnow
        end
        %         set(hLine.MarkerHandle,'EdgeColorBinding','discrete','EdgeColorData',uint8(255*colind)')
        clear hLine
        kk=kk+1;
    end
    
end