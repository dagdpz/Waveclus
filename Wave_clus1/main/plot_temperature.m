function handles=plot_temperature(handles)
cax=handles.axesTemp;
% axes(handles.axesTemp);
my=max(max(handles.tree(:,5:end)))*1.1;
mx=size(handles.tree,1);

semilogy(cax,1:mx,handles.tree(:,5:end));
ylim(cax,[1 my])
xlim(cax,[0.5 mx+0.5]);
handles.hver=line([handles.temp handles.temp],[1 my],'linestyle',':','color','k','parent',cax);
handles.hhor=line([1 mx],[handles.min_clus handles.min_clus],'linestyle',':','color','k','parent',cax);
