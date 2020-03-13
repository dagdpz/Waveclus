function handles=plot_details(handles)

set(handles.mainfig,'Visible','Off');
% set(handles.hdetailsfig,'Visible','On');
figure(handles.hdetailsfig);

handles=get_clustermean(handles);

%temporarily accepted spikes (indices in the array of all spikes, sorted accrding to the distance to the center)
handles.details_inds_rejected=[];%indices of spikes rejected because of other reasons (not the distance to mean)
handles.details_inds_accepted=setdiff(handles.details_inds_ofdistsorted,handles.details_inds_rejected);

t=length(handles.details_inds_accepted);
%initialize sliders
set(handles.hdetailsallspikessliderH,'Min',1);
set(handles.hdetailsallspikessliderH,'Max',t);
set(handles.hdetailsallspikessliderH,'Value',t);
set(handles.hdetailsallspikessliderH,'SliderStep',[1/(t-1) 2/(t-1)]);
handles.hdetailsallspikessliderH_prev=t;

step=1;
set(handles.hdetailsallspikessliderV,'Min',1);
set(handles.hdetailsallspikessliderV,'Max',t);
set(handles.hdetailsallspikessliderV,'Value',step);
% set(handles.hdetailsallspikessliderV,'SliderStep',[ceil(t/100) ceil(t/10)]/t);
set(handles.hdetailsallspikessliderV,'SliderStep',[1/(t-1) ceil(t/10)/t]);
set(handles.hdetailstextStep,'String',num2str(step));

set(handles.hdetailsprune,'Value',0);

handles=details_removedistantspikes(handles);

plot_details_shapeevolution(handles);

plot_details_instFR(handles);
plot_details_distinctfeatures(handles);
plot_details_closebyspikes(handles);

