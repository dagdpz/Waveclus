function handles=wc_plot_spikes_many_clusters(handles,clusterstoplot,ax)
% plots spike shapes for several clusters
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings
%define limits
clusterstoplot(clusterstoplot==0)=numel(handles.classind);
toplot=1;

cla(ax);
hold(ax,'on');
n_ss=length(handles.sp_time);
handles.mean_ss=NaN(max(clusterstoplot),size(handles.spikes,2));
handles.std_ss=NaN(max(clusterstoplot),size(handles.spikes,2));
for i=clusterstoplot;
    if i==length(handles.classind), color= [0 0 0];
    else color=handles.colors(mod(i-1,size(handles.colors,1))+1,:); end
    spikes=handles.spikes(handles.classind{i},:);
    len=size(spikes,1);
    if len==0, cla(ax); title(ax,'');continue; end
    %mean and std for whole cluster independent from plotted number
    clear m s
    m(1,1:n_ss)=mean(spikes,1);
    s(1,1:n_ss)=std(spikes,[],1);
    if isempty(m)
        m=zeros(1,n_ss);
        s=zeros(1,n_ss);
    end
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
    
    if toplot~=2
        if i == size(handles.classind,2)
        else
            plot(ax,handles.sp_time, m,'color',color,'linewidth',2);
            plot(ax,handles.sp_time, m+s,'color',color);
            plot(ax,handles.sp_time, m-s,'color',color);
        end
    else
        plot(ax,handles.sp_time, m,'color',color,'linewidth',2);
    end
end
lims=[min(handles.sp_time) max(handles.sp_time) ...
    min([min(handles.mean_ss(clusterstoplot(2:end),:)-handles.std_ss(clusterstoplot(2:end),:)),0])*1.2 ...
    max([max(handles.mean_ss(clusterstoplot(2:end),:)+handles.std_ss(clusterstoplot(2:end),:)),0.0000001])*1.2];
handles.lims = lims;
axis(ax,lims);
set(ax,'XGrid','on');