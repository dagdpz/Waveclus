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
