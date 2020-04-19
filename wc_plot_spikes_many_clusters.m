function handles=wc_plot_spikes_many_clusters(handles,clusterstoplot,ax)
% plots spike shapes for several clusters
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings
%define limits
clusterstoplot(clusterstoplot==0)=sum(~cellfun(@isempty,handles.classind));
t=handles.spikes([handles.classind{clusterstoplot}],:);
%--------------
% lims=[min(handles.sp_time) max(handles.sp_time) min(min(t)) max(max(t))];
%-----------------
%toplot=get(handles.htoplot,'Value'); %% this needs to go somewhere else
toplot=1;

% axes(ax);
% cla
% hold on
cla(ax);
hold(ax,'on');
n_ss=length(handles.sp_time);
handles.mean_ss=NaN(max(clusterstoplot),size(handles.spikes,2));
handles.std_ss=NaN(max(clusterstoplot),size(handles.spikes,2));
for i=clusterstoplot;
    if i==length(handles.classind), color= [0 0 0];
    else color=handles.colors(i,:); end
    
    spikes=handles.spikes(handles.classind{i},:);

    len=size(spikes,1);
    %     if len==0, cla; title('');continue; end
    if len==0, cla(ax); title(ax,'');continue; end
    %mean and std for whole cluster independent from plotted number
    clear m s 
    m(1,1:n_ss)=mean(spikes,1);
    s(1,1:n_ss)=std(spikes,[],1);
    handles.mean_ss(i,1:n_ss)=m;
    handles.std_ss(i,1:n_ss)=s;

    switch toplot,
        case 1, sp=spikes;%all
        case 2, sp=spikes; %average
        case 3, sp=spikes(round(linspace(1,len,len/20)),:);%  5%
        case 4, sp=spikes(round(linspace(1,len,len/10)),:);% 10%
        case 5, sp=spikes(round(linspace(1,len,len/5)),:);%  20%
        case 6, sp=spikes(round(linspace(1,len,len/2)),:);%  50%
    end
    clear spikes
    
    newlen=size(sp,1);
    tind=1:newlen;
    
    if newlen>MAX_SPIKES_TO_PLOT,
        tind=round(linspace(1,newlen,MAX_SPIKES_TO_PLOT));
        newlen=MAX_SPIKES_TO_PLOT;
    end
    
    sp=sp(tind,:);

    %---------------------
    if toplot~=2
        if i == size(handles.classind,2)
        else
        plot(ax,handles.sp_time, m,'color',color,'linewidth',2);
        plot(ax,handles.sp_time, m+s,'color',color);
        plot(ax,handles.sp_time, m-s,'color',color);
        end
%         set(gcf,'CurrentAxes',ax)
%         fill([handles.sp_time,fliplr(handles.sp_time)],[m+s,fliplr(m-s)],color, 'FaceAlpha', 0.25,'EdgeColor','none');
    
%     if toplot~=2, plot(ax,handles.sp_time, sp','color',color); 
    else plot(ax,handles.sp_time, m,'color',color,'linewidth',2); end
    %     if toplot~=2, plot(handles.sp_time, sp','color',color);
    %     else plot(handles.sp_time, m,'color',color,'linewidth',2); end
    %----------------------
    
end
%-------------------
% v = axis(ax);
lims=[min(handles.sp_time) max(handles.sp_time) min(min(handles.mean_ss(clusterstoplot(2:end),:)-handles.std_ss(clusterstoplot(2:end),:)))*1.2 max(max(handles.mean_ss(clusterstoplot(2:end),:)+handles.std_ss(clusterstoplot(2:end),:)))*1.2];
handles.lims = lims;
axis(ax,lims);
%---------------------
set(ax,'XGrid','on');
% axis(lims);