function handles=wc_create_featuresfig(handles)
if ishandle(handles.mainfig) && isfield(handles,'hfeatures') && isvalid(handles.hfeatures) %gui open, timecoursefigurehandle exists and figure is there
    figure(handles.hfeatures);clf;
else
    handles.hfeatures=figure('Visible','On','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
        'Name','Features','NumberTitle','off','Color','w',...
        'UserData',handles.mainfig,...
        'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
        'Tag','Features');
end

nf=handles.nfeatures;
np=nf;%*(nf-1)/2;
nrow=floor(sqrt(np));
ncol=ceil(sqrt(np));
stepx=0.01; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.01; hight=(1-(nrow+1)*stepy)/nrow;
mnx=prctile(handles.features(:,2),0.1);
mxx=prctile(handles.features(:,2),99.9);
mny=prctile(handles.features(:,2),0.1);
mxy=prctile(handles.features(:,2),99.9);
k=1;
handles.hfaxes(1)=axes('position',[stepx+(width+stepx)*mod(k-1,ncol) 1-(stepy+hight)*(ceil(k/ncol)) width hight],...
    'Tag',sprintf('Cluster%d',k),...
    'xtick',[],'ytick',[],'xlim',[mnx mxx],'ylim',[mny mxy],...
    'Parent',handles.hfeatures,'ButtonDownFcn',{@wc_cluster_cut_callback},'NextPlot','add');
h=title(sprintf('histogram feature: %s',handles.feature_names{2}));
set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
    'Parent',handles.hfaxes(k));
box on
k=2;
for i=2 %:nf,
    for j=i:nf,
        mnx=prctile(handles.features(:,i),0.1);
        mxx=prctile(handles.features(:,i),99.9);
        mny=prctile(handles.features(:,j),0.1);
        mxy=prctile(handles.features(:,j),99.9);
        mnx=mnx-abs(mnx)*0.1;
        mxx=mxx+abs(mxx)*0.1;
        mny=mny-abs(mny)*0.1;
        mxy=mxy+abs(mxy)*0.1;
        handles.hfaxes(k)=axes('position',[stepx+(width+stepx)*mod(k-1,ncol) 1-(stepy+hight)*(ceil(k/ncol)) width hight],...
            'Tag',sprintf('Cluster%d',k),'UserData',handles.mainfig,...
            'xtick',[],'ytick',[],'xlim',[mnx mxx],'ylim',[mny mxy],...
            'Parent',handles.hfeatures,'ButtonDownFcn',{@wc_cluster_cut_callback},'NextPlot','add');
        if j==i
            h=title(sprintf('histogram feature: %s',handles.feature_names{i+1}));
        else
            h=title(sprintf('%s vs %s',handles.feature_names{j},handles.feature_names{i}));
        end
        set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
            'Parent',handles.hfaxes(k));
        box on
        k=k+1;
    end
end

function copy_to_new_window(source,event)
figure(100);
cla;
copyobj(get(source,'children'),gca);