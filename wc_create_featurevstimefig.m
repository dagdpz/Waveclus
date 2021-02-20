function handles=wc_create_featurevstimefig(handles)
if ishghandle(handles.mainfig) && isfield(handles,'htimecourse') && ishghandle(handles.htimecourse) %&& isvalid(handles.htimecourse) %gui open, timecoursefigurehandle exists and figure is there
    figure(handles.htimecourse);clf;
else
    handles.htimecourse=figure('Visible','Off','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
        'Name','Features vs time','NumberTitle','off','Color','w',...
        'UserData',handles.mainfig,...
        'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
        'Tag','Time');
end

nf=handles.nfeatures;
nrow=ceil(sqrt(nf));
ncol=ceil(sqrt(nf));

stepx=0.01; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.01; hight=(1-(nrow+1)*stepy)/nrow;
mny=min(min(handles.index));
mxy=max(max(handles.index));
for i=1:nf,
    mnx=prctile(handles.features(:,i),0.1);
    mxx=prctile(handles.features(:,i),99.9);
    handles.htaxes(i)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*(ceil(i/ncol)) width hight],...
        'Tag',sprintf('Time_vs_feature%d',i-1),'UserData',handles.mainfig,...
        'xtick',[],'ytick',[],'xlim',[mnx*1.1 mxx*1.1],'ylim',[mny mxy],...
        'Parent',handles.htimecourse,'ButtonDownFcn',{@wc_cluster_cut_callback},'NextPlot','add');
    h=title(sprintf('time vs %s',handles.feature_names{i}));
    set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
        'Parent',handles.htaxes(i));
    box on
end

function copy_to_new_window(source,event)
figure(100);
cla;
copyobj(get(source,'children'),gca);