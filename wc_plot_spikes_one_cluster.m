function wc_plot_spikes_one_cluster(handles,clustertoplot)
% plots spike shape for one cluster, cluster zero treated slightly
% differently
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings

%define differences between cluster 0 and other clusters
%define scale on y axes
if clustertoplot==0, %all spikes
    i=length(handles.classind);
    avecolor='c';
    color=[0 0 0];
    ax=handles.spikeaxes(end);
else %except for cluster 0
    i=clustertoplot;
    avecolor='k';
    color=handles.colors(i,:);
    ax=handles.spikeaxes(i);
end

%toplot=get(handles.htoplot,'Value'); %how many spikes to plot
toplot=1; %how many spikes to plot

handles.sp_time=1:(handles.WC.w_pre+handles.WC.w_post);

hold(ax,'on');
spikes=handles.spikes(handles.classind{i},:);
len=size(spikes,1);
if len==0, cla(ax); title(ax,'');return; end

switch toplot,
    case 1, sp=spikes;                                  % all
    case 2, sp=spikes;avecolor=color;                   % average
    case 3, sp=spikes(round(linspace(1,len,len/20)),:); % 5%
    case 4, sp=spikes(round(linspace(1,len,len/10)),:); % 10%
    case 5, sp=spikes(round(linspace(1,len,len/5)),:);  % 20%
    case 6, sp=spikes(round(linspace(1,len,len/2)),:);  % 50%
end

clear spikes
newlen=size(sp,1);
tind=1:newlen;

if newlen>MAX_SPIKES_TO_PLOT,
    tind=round(linspace(1,newlen,MAX_SPIKES_TO_PLOT));
end

sp=sp(tind,:);

if toplot~=2,
    VER=version('-release');
    VER=str2double(VER(1:4));
    if VER>=2014
        plot(ax,handles.sp_time, sp','color',[color 0.1]);
    else
        plot(ax,handles.sp_time, sp','color',color);
    end
end

%mean and std for whole cluster independent from plotted number
m=handles.mean_ss(i,:);
s=handles.std_ss(i,:);

%define width and height of spike, assumes negative threshold
sr=handles.WC.sr;
[w1, w2, amp, ~]=wc_spike_width(m);
w1=w1/sr*1e6;
w2=w2/sr*1e6;

plot(ax,handles.sp_time, m,'color',avecolor,'linewidth',2);
plot(ax,handles.sp_time, m+s,'color',avecolor,'linewidth',0.5);
plot(ax,handles.sp_time, m-s,'color',avecolor,'linewidth',0.5);
if toplot~=2,
    title(ax,sprintf('#%d (%d), w=%1.0f (%1.0f) h=%1.0f',len,size(sp,1),w1,w2,amp));
else
    title(ax,sprintf('#%d (%d), w=%1.0f (%1.0f), h=%1.0f',len,1,w1,w2,amp));
end

lims = handles.lims;
axis(ax,lims);
set(ax,'XGrid','on');
