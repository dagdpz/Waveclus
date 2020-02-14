function plot_spikes(ax,time,sp,lims,color,toplot,max_spikes_to_plot,superimposed)
if exist('superimposed','var'), sup=1; else sup=0; end
if ~exist('max_spikes_to_plot','var'), max_spikes_to_plot=1000; end

hold(ax,'on');
% axes(ax);hold on;
len=size(sp,1);
if color=='k', c='c'; else c='k'; end
switch toplot,
    case 2, c=color;
    case 3, sp=sp(round(linspace(1,len,len/20)),:);
    case 4, sp=sp(round(linspace(1,len,len/10)),:);
    case 5, sp=sp(round(linspace(1,len,len/5)),:);
    case 6, sp=sp(round(linspace(1,len,len/2)),:);
end
if toplot~=2, plot(ax,time, sp(1:min(end,max_spikes_to_plot),:),'color',color); end
m=mean(sp);s=std(sp);
if (~sup), 
    plot(ax,time, m,'color',c,'linewidth',2);
    plot(ax,time, m+s,'color',c,'linewidth',0.5);
    plot(ax,time, m-s,'color',c,'linewidth',0.5);
    title(ax,sprintf('#%d',size(sp,1)));
elseif (toplot==2)
    plot(ax,time, m,'color',c,'linewidth',2);
end
axis(ax,lims);
set(ax,'XGrid','on');
