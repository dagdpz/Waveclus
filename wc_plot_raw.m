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
