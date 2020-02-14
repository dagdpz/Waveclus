function handles=details_removedistantspikes(handles)
%changes the value of horizontal slider (number os spikes in cluster)
%this procedure is called when 
%handles.hdetailstoplot
%handles.hdetailsallspikessliderH
%handles.hdetailsallspikessliderV
%are changed



[t,ind]=setdiff(handles.details_inds_ofdistsorted,handles.details_inds_rejected);
%setdiff sorts output but we want to save an order accoring to the distance to the mean
handles.details_inds_accepted=[handles.details_inds_ofdistsorted(sort(ind))];


step=round(get(handles.hdetailsallspikessliderV,'Value'));
prevvalue=handles.hdetailsallspikessliderH_prev;

handles.details_inds_accepted=handles.details_inds_accepted(1:prevvalue);

newvalue=get(handles.hdetailsallspikessliderH,'Value');
if abs(newvalue-prevvalue)==2,%big step, until next closeby
    %find most distant closeby spike
    %maybe take 
    MINISI=handles.MINISI;
    MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings
    tt=get(handles.hdetailsclosebyspikesslider,'Value')+1;
    tt=min(tt,length(MINISI));
    minisi=MINISI(tt);

    
    inds=handles.details_inds_accepted;
    inds=sort(inds);%sort indices in time
    timeinds=handles.index(inds);%size is 1:newvalue, values are times of spikes of closebyspikes
    d=diff(timeinds);
    dd=find(d<=minisi);
    %check whether there are any closeby spikes, otherwise proceed as for small step 
    if isempty(dd), newvalue=prevvalue+sign(newvalue-prevvalue)*step;
    else
        [tt,ia,ib]=intersect(handles.details_inds_accepted,inds([dd;dd+1]));        
        newvalue=max(ia)-1;
    end
else %small step, step based on vertical slider
    newvalue=prevvalue+sign(newvalue-prevvalue)*step;
end

maxvalue=length(t);
if newvalue>maxvalue, newvalue=maxvalue; end
minvalue=get(handles.hdetailsallspikessliderH,'Min');
if newvalue<minvalue, newvalue=minvalue; end

set(handles.hdetailsallspikessliderH,'Value',newvalue);
handles.hdetailsallspikessliderH_prev=newvalue;



handles.details_inds_accepted=[handles.details_inds_ofdistsorted(sort(ind))];
handles.details_inds_accepted=handles.details_inds_accepted(1:newvalue);

%update handles
% Update handles structure
guidata(handles.mainfig, handles);

plot_details_distantspikes(handles);