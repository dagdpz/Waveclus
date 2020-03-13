function plot_details_closebyspikes(handles)
%procedure which plots close by spikes

MINISI=handles.MINISI;

MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings
t=get(handles.hdetailsclosebyspikesslider,'Value')+1;
t=min(t,length(MINISI));
minisi=MINISI(t);

inds=handles.details_inds_accepted;
%size is 1:newvalue, values are indices in array of all indices
%indices are sorted according to the distance to the mean

inds=sort(inds);%sort indices in time

timeinds=handles.index(inds);%size is 1:newvalue, values are times of spikes of closebyspikes

d=diff(timeinds);
dd=find(d<=minisi);
%values are indices of spike times of closeby spikes in array of current cluster
% axes(handles.hdetailsclosebyspikes); cla; hold on;
cax=handles.hdetailsclosebyspikes;
cla(cax);
hold(cax,'on');

cinds=inds(dd);
if ~isempty(cinds), 
    %     set(handles.hdetailsclosebyspikes,'Visible','on');
    %     plot(handles.sp_time,handles.spikes(unique(cinds(round(linspace(1,length(cinds),MAX_SPIKES_TO_PLOT)))),:)','b');
    set(cax,'Visible','on');
    plot(cax,handles.sp_time,handles.spikes(unique(cinds(round(linspace(1,length(cinds),MAX_SPIKES_TO_PLOT)))),:)','b');
end
cinds=inds(dd+1);
if ~isempty(cinds), 
    set(cax,'Visible','on');
    plot(cax,handles.sp_time,handles.spikes(unique(cinds(round(linspace(1,length(cinds),MAX_SPIKES_TO_PLOT)))),:)','r');
%     set(handles.hdetailsclosebyspikes,'Visible','on');
%     plot(handles.sp_time,handles.spikes(unique(cinds(round(linspace(1,length(cinds),MAX_SPIKES_TO_PLOT)))),:)','r');
end

xlim(cax,[min(handles.sp_time) max(handles.sp_time)]);
t=get(cax,'ylim');
text(min(handles.sp_time),t(2),sprintf('%d<%1.2fms',length(dd),minisi),...
    'horizontalalignment','left','verticalalignment','bottom','parent',cax);


% handles.currdist(unique([dd dd+1])) gives distances of closeby spikes 