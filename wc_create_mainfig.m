function handles=wc_create_mainfig(handles)
handles.mainfig=figure('Visible','on','Units','Normalized','Position',[0 0 1,0.9],...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Name','Wave_clus','NumberTitle','off',...
    'Tag','Mainfig');

figure(handles.mainfig);
ncl=handles.ncl;
ncol=5;
nrow=3;
if ncl<(ncol-1)*2
nrow=2;
end

stepx=0.02; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;

handles.axesAllClusters=axes('position',[stepx 1-(stepy+hight)*( 1 ) width hight],...
    'Tag','AllClusters','NextPlot','add');
handles.axesClust0=axes('position',[stepx+(width+stepx)*4 1-(stepy+hight)*( nrow ) width hight],...
    'Tag','Clust0','NextPlot','add');
handles.axesTemp=axes('position',[stepx 1-(stepy+hight)*( 2 )+0.1*hight width hight],...
        'Tag','Temperature');
    
handles.spikeaxes=[];
for i=1:nrow*(ncol-1)-1
    r=ceil(i/(ncol-1));
    c=mod(i-1,ncol-1)+1;
    handles.spikeaxes(i)=axes('position',[stepx+(width+stepx)*c 1-(stepy+hight)*( r ) width hight],...
        'Tag',sprintf('Clust%d',i),'NextPlot','add');
end
handles.spikeaxes(end+1)=handles.axesClust0;
