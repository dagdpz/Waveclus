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