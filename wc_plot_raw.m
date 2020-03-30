function wc_plot_raw(handles)
%add threshold
cax=handles.axesTS;
cla(cax);
line(handles.ts_time, handles.ts,'color','k','parent',cax);
line(handles.ts_time, repmat(-handles.WC.thr,[length(handles.ts_time) 1]),'color','r','parent',cax);
line(handles.ts_time, repmat(handles.WC.thr,[length(handles.ts_time) 1]),'color','r','parent',cax);
set(cax,'xlim',[min(handles.ts_time) max(handles.ts_time)]);
set(cax,'ylim',[min(handles.ts) max(handles.ts)]*1.05);
hold(cax,'on');

%------
% for i=1:handles.ncl,
%     t=find(handles.index(handles.classind{i})<max(handles.ts_time)*1000);
%     plot(cax,handles.index(handles.classind{i}(t))/1000,double(max(handles.ts))*ones(1,length(t)),'color',handles.colors(i),'linestyle','none','marker','.','markersize',5);
% end
%------


% axes(handles.axesTS);
% cla
% line(handles.ts_time, handles.ts,'color','k');
% set(gca,'xlim',[min(handles.ts_time) max(handles.ts_time)]);
% set(gca,'ylim',[min(handles.ts) max(handles.ts)]*1.05);
% hold on
% for i=1:handles.ncl,
%     t=find(handles.index(handles.classind{i})<max(handles.ts_time)*1000);
%     
%     plot(handles.index(handles.classind{i}(t))/1000,max(handles.ts)*ones(1,length(t)),'color',handles.colors(i),'linestyle','none','marker','.','markersize',5);
% end
