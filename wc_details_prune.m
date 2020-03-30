function handles=wc_details_prune(handles)

newvalue=get(handles.hdetailsallspikessliderH,'Value');
if get(handles.hdetailsprune,'Value')==1,
    
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
    cinds=union(inds(dd),inds(dd+1));%values are indices of closeby spikes in array of all indices 
    
    %reject based on the distance to the cluster center
    reject_inds=[];
    for i=1:length(dd),
        if handles.details_dist2center(dd(i))>handles.details_dist2center(dd(i)+1),
            reject_inds=[reject_inds inds(dd(i))];
        else
            reject_inds=[reject_inds inds(dd(i)+1)];
        end        
    end
    handles.details_inds_rejected=reject_inds;
    newvalue=newvalue-length(handles.details_inds_rejected);
    
    %%%%plot rejected
    figure(100);
    clf; hold on;
    %this is for legend
    line([0],[0],'color','r','Visible','off');
    line([0],[0],'color','b','Visible','off');
    
    plot(handles.sp_time,handles.spikes(reject_inds,:),'r');
    hold on
    plot(handles.sp_time,handles.spikes(setdiff(cinds,reject_inds),:),'b');
    legend({'Rejected','Left'})
    
else
    newvalue=newvalue+length(handles.details_inds_rejected);
    handles.details_inds_rejected=[];
end

handles.hdetailsallspikessliderH_prev=newvalue;
set(handles.hdetailsallspikessliderH,'Value',newvalue);

handles=details_removedistantspikes(handles);
plot_details_shapeevolution(handles);
plot_details_instFR(handles);
plot_details_distinctfeatures(handles);
plot_details_closebyspikes(handles);

% handles.currdist(unique([dd dd+1])) gives distances of closeby spikes 