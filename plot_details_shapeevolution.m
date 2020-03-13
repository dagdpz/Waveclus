function plot_details_shapeevolution(handles)

inds=handles.details_inds_accepted;
inds=sort(inds);

nspk=handles.nspk;
mn=min(min(handles.spikes(inds,:)));
mx=max(max(handles.spikes(inds,:)));
lims=[min(handles.sp_time) max(handles.sp_time) mn mx];
%plot how spike shapes change over time
%all spikes
cla(handles.hdetailsaveragespikesaxes);hold(handles.hdetailsaveragespikesaxes,'on');mx=-Inf;mn=Inf;
cla(handles.hdetailsisiaxes);hold(handles.hdetailsisiaxes,'on');
% inds=handles.classind{handles.details_currcluster};

timeinds=handles.index(inds);
time0=min(handles.index);
timeend=max(handles.index);
time5=linspace(time0,timeend,6);

toplot=get(handles.hdetailstoplot,'Value');
MAX_SPIKES_TO_PLOT=handles.const_MAX_SPIKES_TO_PLOT; %to prevent large plottings
for i=1:5,
    inds5=find( (timeinds>=time5(i)) & (timeinds<time5(i+1)));
    cax=handles.hdetailsspikeaxes(i);
    cla(cax);
    if isempty(inds5), title('');continue; end
    %     axes(handles.hdetailsspikeaxes(i));cla;
    spikes=handles.spikes(inds(inds5),:);
    
    [w1, w2, amp, halfamp]=spike_width(mean(spikes,1));
    w1=w1/handles.par.sr*1e6;
    w2=w2/handles.par.sr*1e6;

    
    len=size(spikes,1);
    if ~len, title(sprintf('#0 (0)')); continue; end
    switch toplot,
        case 1, sp=spikes;%all
        case 2, sp=spikes(round(linspace(1,len,max(1,len/20))),:);
        case 3, sp=spikes(round(linspace(1,len,max(1,len/10))),:);
        case 4, sp=spikes(round(linspace(1,len,max(1,len/5))),:);
        case 5, sp=spikes(round(linspace(1,len,max(1,len/2))),:);
    end

    newlen=size(sp,1);
    tind=1:newlen;

    if newlen>MAX_SPIKES_TO_PLOT,
        tind=round(linspace(1,newlen,MAX_SPIKES_TO_PLOT));
        newlen=MAX_SPIKES_TO_PLOT;
    end
    
    plot(cax,handles.sp_time,sp(tind,:)','color',handles.colors(i));
    %     plot(handles.sp_time,sp','color',handles.colors(i));
    %     axis(lims);
    axis(cax,lims);
    title(cax,sprintf('#%d (%d), w=%1.0f (%1.0f) h=%1.0f',len,newlen,w1,w2,amp));
    %     title(cax,sprintf('#%d (%d)',len,newlen),'fontsize',12);
    %isi
    cax=handles.hdetailsisiaxes;
    edges=0:1:100;
    edges1=0:1:100;
    if length(inds(inds5))>1,
        [N,X]=histc(diff(handles.index(inds(inds5))),edges);
        %         if ~all(N==0), plot(edges1,spline(edges,N/max(N),edges1),'color',handles.colors(i)); end
        if ~all(N==0), plot(cax,edges,smooth(N/max(N),5),'color',handles.colors(i)); end
    end
    %     %isi
    %     axes(handles.hdetailsisiaxes);hold on
    %     edges=0:1:100;
    %     edges1=0:1:100;
    %     if length(inds(inds5))>1,
    %         [N,X]=histc(diff(handles.index(inds(inds5))),edges);
    %         %         if ~all(N==0), plot(edges1,spline(edges,N/max(N),edges1),'color',handles.colors(i)); end
    %         if ~all(N==0), plot(edges,smooth(N/max(N),5),'color',handles.colors(i)); end
    %     end
    
    %averages
    cax=handles.hdetailsaveragespikesaxes;
    mm=mean(handles.spikes(inds(inds5),:),1);
    plot(cax,handles.sp_time,mm,'color',handles.colors(i));
    t=min(mm); if t<mn, mn=t; end
    t=max(mm); if t>mx, mx=t; end
    %     %averages
    %     axes(handles.hdetailsaveragespikesaxes);hold on
    %     mm=mean(handles.spikes(inds(inds5),:),1);
    %     plot(handles.sp_time,mm,'color',handles.colors(i));
    %     t=min(mm); if t<mn, mn=t; end
    %     t=max(mm); if t>mx, mx=t; end
end
cax=handles.hdetailsaveragespikesaxes;
% axes(handles.hdetailsaveragespikesaxes);hold on
ylim(cax,[mn mx]*1.05);
xlim(cax,[min(handles.sp_time) max(handles.sp_time)]);


