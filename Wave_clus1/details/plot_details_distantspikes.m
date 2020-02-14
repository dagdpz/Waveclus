function plot_details_distantspikes(handles)
%plot the most distant from cenetr spikes
%highlight in red the most most distant (candidates to delete)
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings

step=round(get(handles.hdetailsallspikessliderV,'Value'));%number of red spikes, first candidates to be removed
inds=handles.details_inds_accepted;

spikes=handles.spikes(inds,:);
lims=[min(handles.sp_time) max(handles.sp_time) min(spikes(:)) max(spikes(:))];

toplot=get(handles.hdetailstoplot,'Value');

% axes(handles.hdetailsallspikes);
% cla; hold on
cax=handles.hdetailsallspikes;
cla(cax); hold(cax,'on');

%how many spikes to plot
len=size(spikes,1);
switch toplot,
    case 1, sp=spikes;%all
    case 2, sp=spikes(end-round(len/20):end,:);%%5
    case 3, sp=spikes(end-round(len/10):end,:);%%10
    case 4, sp=spikes(end-round(len/5):end,:);%%20
    case 5, sp=spikes(end-round(len/2):end,:);%%50
end

sp=sp(end-min(end,MAX_SPIKES_TO_PLOT)+1:end,:);

%plot all spikes (not exactly all, defined by toplot)
plot(cax,handles.sp_time, sp','color','b');

plot(cax,handles.sp_time, handles.details_mean_spikesshape,'color','k','linewidth',2);
plot(cax,handles.sp_time,sp(max(end-step+1,1):end,:)','r','linewidth',1);
title(cax,sprintf('#%d (%d) (%d)',size(sp,1),len,length(handles.details_dist2center)));
axis(cax,lims);