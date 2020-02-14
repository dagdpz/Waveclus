function handles=plot_spikes_isi(handles)
ncl=handles.ncl;

%clear plotted, do not clear all clusters after rejection operation only 
plotted=find(handles.plotted);
for i=plotted,
    if i>=handles.rejected,
        cla(handles.spikeaxes(i));
        cla(handles.isiaxes(i));
        %         axes(handles.spikeaxes(i));cla;
        %         axes(handles.isiaxes(i));cla;
        set(handles.hclustergroup{i},'Visible','off');%make axes and assoiciated buttons invisible
    end
end
%cluster zero
cla(handles.spikeaxes(end));
cla(handles.isiaxes(end));
% axes(handles.spikeaxes(end));cla;
% axes(handles.isiaxes(end));cla;

if handles.nspk~=size(handles.spikes,1), error('Number of spikes is different in spikes and in spc output file'); end
%all spikes superimposed
handles=plot_spikes_many(handles,0:ncl,handles.axesAllClusters);


for i=handles.rejected:ncl,
    set(handles.hclustergroup{i},'Visible','on');
    plot_spikes_one(handles,i);
    plot_isi(handles, i);
end

%cluster zero
plot_spikes_one(handles,0);
plot_isi(handles,0);

% if length(handles.classind{end}), 
%     mn=min(min(handles.spikes(handles.classind{end},:)));
%     mx=max(max(handles.spikes(handles.classind{end},:)));
%     lims=[min(handles.sp_time) max(handles.sp_time) mn mx];
%     if length(handles.classind{end}), 
%         plot_spikes(handles.spikeaxes(end),...
%             handles.sp_time,handles.spikes(handles.classind{end},:),...
%             lims,'k',get(handles.htoplot,'Value'));
%         plot_isi(handles.isiaxes(end),...
%             diff(handles.index(handles.classind{end})),...
%             'k',get(handles.htoplot,'Value'));
%     end
% else
%     cla(handles.spikeaxes(end));
%     cla(handles.isiaxes(end));
%     title('');
% end

if ncl>3, set(handles.hsuppl,'Visible','on'); 
else set(handles.hsuppl,'Visible','off'); 
end

handles.plotted=[];
handles.plotted(1:ncl)=1;
axes(handles.axesTemp)