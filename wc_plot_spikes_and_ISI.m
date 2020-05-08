function handles=wc_plot_spikes_and_ISI(handles)
ncl=handles.ncl;

%clear plotted, do not clear all clusters after rejection operation only 
plotted=find(handles.plotted);
for i=plotted,
    if i>=handles.rejected,
        cla(handles.spikeaxes(i));
        if handles.isaGUI
            set(handles.hclustergroup{i},'Visible','off');%make axes and assoiciated buttons invisible
        end
    end
end
%cluster zero
cla(handles.spikeaxes(end));
if handles.nspk~=size(handles.spikes,1), disp('Number of spikes is different in spikes and in spc output file'); end
%all spikes superimposed

handles.sp_time=1:(handles.WC.w_pre+handles.WC.w_post);
handles=wc_plot_spikes_many_clusters(handles,0:ncl,handles.axesAllClusters);

for i=handles.rejected:ncl,
    if handles.isaGUI
        set(handles.hclustergroup{i},'Visible','on');
    end
    wc_plot_spikes_one_cluster(handles,i);
    wc_plot_ISI(handles, i);
end

%cluster zero
wc_plot_spikes_one_cluster(handles,0);
wc_plot_ISI(handles,0);

handles.plotted=[];
handles.plotted(1:ncl)=1;

function wc_plot_spikes_one_cluster(handles,clustertoplot)
% plots spike shape for one cluster, cluster zero treated slightly
% differently
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings

%define differences between cluster 0 and other clusters
%define scale on y axes
if clustertoplot==0, %all spikes
    i=sum(~cellfun(@isempty,handles.classind));
    avecolor=[0.5 0.5 0.5];
    color=[0 0 0];
    ax=handles.spikeaxes(end);
else %except for cluster 0
    i=clustertoplot;
    avecolor='k';
    color=handles.colors(mod(i-1,size(handles.colors,1))+1,:);
    ax=handles.spikeaxes(i);
end

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

function wc_plot_ISI(handles,i)
if i==0, 
    i = length(handles.classind); 
    ax=handles.spikeaxes(end); 
    color='k';
else
    ax=handles.spikeaxes(i);     
    color=handles.colors(mod(i-1,size(handles.colors,1))+1,:);
end
isi=diff(handles.index(handles.classind{i}));

%cla(ax);
edges=0:1:100;
if isempty(isi), isi=0; end
N=histc(isi,edges);
x_lim=get(ax,'xlim');
y_lim=get(ax,'ylim');

%% background
x1=x_lim(1)+diff(x_lim)/2;
y1=y_lim(1);
w=diff(x_lim)/2;
h=diff(y_lim)/2;
ca=gca;
axes(ax);
R=rectangle('Position', [x1 y1 w h]);
text(x1+diff(x_lim)/20,y_lim(2)-diff(y_lim)*11/20,sprintf('%d in <2ms, %d in <1ms',sum(N(1:2)),N(1)));
axes(ca); 

set(R,'facecolor','w');
N_transformed=N/max(N)*diff(y_lim)/2;%+y_lim(1);
edges_transformed=edges/max(edges)*diff(x_lim)/2+diff(x_lim)/2+x_lim(1);
h=bar(ax,edges_transformed,N_transformed,'histc');
ytoshift=get(h,'ydata');
ytoshift=ytoshift+y_lim(1);
set(h,'ydata',ytoshift);
set(h,'facecolor',color,'edgecolor',color,'linewidth',0.01);    