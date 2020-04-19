function wc_plot_ISI(handles,i)

if i==0, 
    i = length(handles.classind); ax=handles.isiaxes(end); color='k';
else
    ax=handles.isiaxes(i); color=handles.colors(i,:);
end
isi=diff(handles.index(handles.classind{i}));

cla(ax);
edges=0:1:100;
if isempty(isi), isi=0; end
[N,X]=histc(isi,edges);
h=bar(ax,edges,N,'histc');
set(h,'facecolor',color,'edgecolor',color,'linewidth',0.01);    
xlim(ax,[0 100]);
title(ax,sprintf('%d in <2ms, %d in <1ms',sum(N(1:2)),N(1)));