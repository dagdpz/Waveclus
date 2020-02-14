function plot_spikes_one(handles,clustertoplot)
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
    t=handles.spikes([handles.classind{1:end}],:);
else %except for cluster 0
    i=clustertoplot;
    avecolor='k';
    color=handles.colors(i,:);
    ax=handles.spikeaxes(i);
    t=handles.spikes([handles.classind{1:end-1}],:);
end

%--------------------------------------
% lims=[min(handles.sp_time) max(handles.sp_time) min(min(t)) max(max(t))];
lims = handles.lims;
%-------------------------------------
%how many spikes to plot
toplot=get(handles.htoplot,'Value');

% axes(ax);
% hold on
hold(ax,'on');
spikes=handles.spikes(handles.classind{i},:);
len=size(spikes,1);
if len==0, cla(ax); title(ax,'');return; end

switch toplot,
    case 1, sp=spikes;%all
    case 2, sp=spikes;avecolor=color; %average
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
%if toplot~=2, plot(ax,handles.sp_time, sp','color',color); end
if toplot~=2, plot(ax,handles.sp_time, sp','color',[color 0.1]); end
%mean and std fro whole cluster independent from plotted number
m=handles.mean_ss(i,:);
s=handles.std_ss(i,:);

%define width and height of spike, assumes negative threshold
w_pre=handles.par.w_pre;
bl=round(w_pre/2);
sr=handles.par.sr;
nstd=3;

[w1, w2, amp, halfamp]=spike_width(m);
w1=w1/sr*1e6;
w2=w2/sr*1e6;

% [mx1,indmx1]=max(m(1:w_pre));
% [mx2,indmx2]=max(m(w_pre+1:end));
% indmx2=indmx2+w_pre;
% 
% t1=m(1:bl);
% st=std(t1); me=mean(t1);
% if (indmx1<bl)||(mx1<me+nstd*st),
% %     if indmx1<bl, warning('w:EarlyMax',...
% %             'Unit %s, ch %s, first maximum is quite early, probably there is none',legclust{j},legfilt{i});
% %     end
% %     if mx1<me+nstd*st, warning('w:SmallMax',...
% %             'Unit %s, ch %s, first maximum is quite small, probably there is none',legclust{j},legfilt{i});
% %     end
%     t2=find(m(1:w_pre)<(me-nstd*st));
%     if isempty(t2), warning('w:w','Strange'); indmx1=1;
%     else indmx1=t2(1); end
% end
% 
% width=round((indmx2-indmx1)/sr*1e6);
% w1=round((indmx2-w_pre)/sr*1e6);
% height=round(max(m)-min(m));


plot(ax,handles.sp_time, m,'color',avecolor,'linewidth',2);
plot(ax,handles.sp_time, m+s,'color',avecolor,'linewidth',0.5);
plot(ax,handles.sp_time, m-s,'color',avecolor,'linewidth',0.5);
if toplot~=2,
    title(ax,sprintf('#%d (%d), w=%1.0f (%1.0f) h=%1.0f',len,size(sp,1),w1,w2,amp));
else
    title(ax,sprintf('#%d (%d), w=%1.0f (%1.0f), h=%1.0f',len,1,w1,w2,amp));
end
% line([handles.sp_time(w_pre) handles.sp_time(w_pre)],lims(3:4),'color','k');
% line([handles.sp_time(indmx2) handles.sp_time(indmx2)],lims(3:4),'color','k');
% line([handles.sp_time(indmx1) handles.sp_time(indmx1)],lims(3:4),'color','k');
axis(ax,lims);
set(ax,'XGrid','on');
