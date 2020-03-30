function handles=wc_plot_details(handles)

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

function handles=get_clustermean(handles)
%estimates a cluster mean in predifined proximity space
%also estimates the distances from each spike in a cluster to the mean

i=handles.details_currcluster;
%do some preprocessing for current cluster
%find center of cluster in feature or spike shape space
if handles.forced(i), inds=handles.classind_unforced{i};
else inds=handles.classind{i}; end

handles.details_mean_spikesshape=mean(handles.spikes(inds,:),1);


switch handles.classify_space,
    case 'features',%all features weighted equally
        %if cluster was forced (aka template matching was used) then use an
        %unforced version to estimate the center
        if handles.forced(i), inds=handles.classind_unforced{i};
        else inds=handles.classind{i}; end

        %estimate center
        mt=mean(handles.features(inds,:),1);
        %estimate distance to the center
        inds=handles.classind{i};
        t=handles.features(inds,:);
        st=[length(inds) 1];
        t=(t-repmat(mt,st))./repmat(std(t,[],1),st);
        handles.details_dist2center=t;
    case 'selectedfeatures',%NOT DONE%weight features according to their importance for current cluster
        %find the most separating features
        if handles.forced(i), inds=handles.classind_unforced{i};
        else inds=handles.classind{i}; end
        inds0=handles.classind{end};
        [nspk,nf]=size(handles.features);
        notinds=setdiff(1:nspk,[inds inds0]);
        if length(notinds)>1,
            for j=1:nf,
                zval(j)=NaN;
                [prs(j),h,stats]=ranksum(handles.features(inds,j),handles.features(notinds,j));
                if isfield(stats,'zval'), zval(j)=-abs(stats.zval); end
            end
            [t,indrs]=sort(zval);
        end
        %normalize featues
        t=handles.features(handles.classind{i},:);
        st=[size(t,1) 1];
        mt=mean(t,1);
        t=(t-repmat(mt,st))./repmat(std(t,[],1),st);
        handles.details_dist2center=t;
    case 'spikeshapes',
        %take whole spike shape
        %         mt=mean(handles.spikes(handles.classind{i},:),1);
        %         handles.details_dist2center=handles.spikes(handles.classind{i},:)-repmat(mt,[length(handles.classind{i}) 1]);
        
        %take only half around min(max)
        ints=handles.par.w_pre-handles.par.w_pre/2:handles.par.w_pre+handles.par.w_post/2-1;
        %         mt=mean(handles.spikes(handles.classind{i},ints),1);
        mt=mean(handles.spikes(handles.classind_unforced{i},ints),1);
        handles.details_dist2center=handles.spikes(handles.classind{i},ints)-repmat(mt,[length(handles.classind{i}) 1]);

end
%square distance to the center of cluster
handles.details_dist2center=sum(handles.details_dist2center.^2,2);
[t,ind]=sort(handles.details_dist2center);
%indices of spikes sorted according to the distance to the center of cluster 
handles.details_inds_ofdistsorted=handles.classind{i}(ind);


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
    
    [w1, w2, amp, halfamp]=wc_spike_width(mean(spikes,1));
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

function plot_details_instFR(handles)

inds=handles.details_inds_accepted;
inds=sort(inds);

nspk=handles.nspk;

timeinds=handles.index(inds);
time0=min(handles.index);
timeend=max(handles.index);
time5=linspace(time0,timeend,6);


cax=handles.hdetailshistogram;
cla(cax); hold(cax,'on')
% axes(handles.hdetailshistogram);
% cla(handles.hdetailshistogram);hold on;

ind=handles.index(inds);
edges=linspace(time5(1),time5(end),300)/1000;%in seconds
[N,X]=histc(ind/1000,edges);
N=N/(edges(2)-edges(1));
plot(cax,edges,N,'color',handles.colors(2));%spkies per second

xlim(cax,[min(handles.index) max(handles.index)]/1000);
ylim(cax,[0 max(N)*1.10]);
title(cax,sprintf('%s, class %d',handles.filename,handles.details_currcluster),'fontsize',12,'interpreter','none');
%%plotting of some kinf of events
% bname=sprintf('%s.mat',handles.bname);
% if exist(bname,'file'), 
%     q=load(bname);
%     len=[0 q.len];
%     for i=1:length(q.ids),
%         line([len(i+1) len(i+1)],[0 max(N)*1.10],'linewidth',2);
%         text((len(i+1)+len(i))/2,max(N),q.strids{i},...
%             'fontsize',16,'verticalalignment','bottom','horizontalalignment','center');
%     end
% end
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

function plot_details_distinctfeatures(handles)
nspk=handles.nspk;
inds=handles.details_inds_accepted;
inds=sort(inds);
timeinds=handles.index(inds);
time0=min(handles.index);
timeend=max(handles.index);
time5=linspace(time0,timeend,6);

%find the most separating features
inds0=handles.classind{end};
nf=size(handles.features,2);
notinds=setdiff(1:nspk,[inds inds0]);
indrs=[1:nf];
zval(1:nf)=NaN;
prs(1:nf)=NaN;
if length(notinds)>1,
    for j=1:nf,
        [prs(j),h,stats]=ranksum(handles.features(inds,j),handles.features(notinds,j));
        if isfield(stats,'zval'), zval(j)=-abs(stats.zval); end
    end
    [t,indrs]=sort(zval);
end
%plot features
k=1;
for i=1:3,
    cax=handles.hdetailsfeatures(k);
    cla(cax);hold(cax,'on');
    %     axes(handles.hdetailsfeatures(k)); cla; hold on;
    mnx=+Inf; mxx=-Inf;
    mny=+Inf; mxy=-Inf;
    x=handles.features(:,indrs(2*i-1));  
    y=handles.features(:,indrs(2*i));
    if ~isempty(inds0), plot(cax,x(inds0),y(inds0),'linestyle','none','markersize',.5,'color','k'); end
    for ii=1:5,
        inds5=find( (timeinds>=time5(ii)) & (timeinds<time5(ii+1)));
        if length(inds5)<1, continue; end
        xx=x(inds(inds5));yy=y(inds(inds5));
        mnx=min(mnx,min(xx)); mxx=max(mxx,max(xx));
        mny=min(mny,min(yy)); mxy=max(mxy,max(yy));
        if ~isempty(inds), plot(cax,xx,yy,'linestyle','none','marker','.','markersize',.5,'color',handles.colors(ii)); end
        h=simple_ellipse(cax,mean(xx),mean(yy),3*std(xx),3*std(yy));
        set(h,'linestyle','-','color',handles.colors(ii),'linewidth',2);
        set(cax,'xtick',[],'ytick',[]);
        xlabel(cax,sprintf('p<%1.3f, z=%1.1f',prs(indrs(2*i-1)),-zval(indrs(2*i-1))));
        ylabel(cax,sprintf('p<%1.3f, z=%1.1f',prs(indrs(2*i)),-zval(indrs(2*i))));
    end            
    if ~isempty(notinds), plot(cax,x(notinds),y(notinds),'linestyle','none','marker','.','markersize',.5,'color','y'); end
    
    axis(cax,[mnx mxx mny mxy]);
    h=title(cax,sprintf('%s vs %s',handles.feature_names{indrs(2*i-1)},handles.feature_names{indrs(2*i)}));
    set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left')
    k=k+1;
end





function h = simple_ellipse(cax,x,y,rx,ry) %simple means that ellipse axes are parallel to x and y 
th = 0:pi/50:2*pi;
xunit = rx * cos(th) + x;
yunit = ry * sin(th) + y;
h = plot(cax,xunit, yunit);
% handles.currdist(unique([dd dd+1])) gives distances of closeby spikes 
