function handles=wc_plot_temperature(handles)
%temperature plot
cax=handles.axesTemp;
subplot(cax);

if isfield(handles,'hcol')
    for i=1:numel(handles.hcol)
        delete(handles.hcol(i));
    end
end
if isfield(handles,'hver')
    for i=1:numel(handles.hver)
        if isobject(handles.hver(i))
        delete(handles.hver(i));
        end
    end
end
handles.hver=[];
handles.hcol=[];

hold off;
tree=handles.tree;
min_clus=handles.min_clus;
%colormap(handles.colors(size(tree,2)-4,:))
semilogy(1:handles.WC.num_temp,tree(2:handles.WC.num_temp+1,5:end));
hold on;
colidx=1;
for c=1:size(handles.WC.clus_per_temp,2)
    x=handles.WC.clus_per_temp(1,c)-1;
    handles.hver(c)=line([x x],[1 max(max(tree(:,5:end)))*1.1],'linestyle',':','color','k');
    %n_clus_t=numel(handles.WC.n_clus_per_temp{c});
    %for c=1:n_clus_t
        cc=handles.WC.clus_per_temp(2,c);
        y=tree(x+1,cc+4);
        handles.hcol(colidx)=scatter(x,y,50,handles.colors(mod(colidx-1,size(handles.colors,1))+1,:),'o','filled');
        colidx=colidx+1;
    %end
end
handles.hhor=line([1 handles.WC.num_temp],[min_clus min_clus],'linestyle',':','color','k');
ylim([1 max(max(tree(:,5:end)))*1.1])
xlim([0.5 handles.WC.num_temp+0.5]);
