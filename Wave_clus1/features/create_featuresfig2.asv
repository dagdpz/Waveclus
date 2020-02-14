function handles=create_featuresfig2(handles);

if isfield(handles,'hfeatures'),
    figure(handles.hfeatures);clf;
else
    if ishandle(handles.mainfig)
        
        handles.hfeatures=figure('Visible','Off','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
            'Name','Features','NumberTitle','off','Color','w',...
            'UserData',handles.mainfig,...
            'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
            'CloseRequestFcn',{@closefeatures_callback},...
            'Tag','Features');
    else %no gui opend
        
        handles.hfeatures=figure('Visible','Off','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
            'Name','Features','NumberTitle','off','Color','w',...
            'UserData',handles.mainfig,...
            'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
            'Tag','Features');
    end
    %     'Paperunits','Normalized','Paperorientation','Landscape','PaperPosition',[0.01 0.01 0.98 0.98],...
    
end

nf=handles.nfeatures;
np=nf*(nf-1)/2;
nrow=floor(sqrt(np));
ncol=ceil(sqrt(np));

stepx=0.01; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.01; hight=(1-(nrow+1)*stepy)/nrow;
k=1;
for i=1:nf,
    for j=i+1:nf,
        mnx=min(min(handles.features(:,i)));
        mny=min(min(handles.features(:,j)));
        mxx=max(max(handles.features(:,i)));
        mxy=max(max(handles.features(:,j)));
        handles.hfaxes(k)=axes('position',[stepx+(width+stepx)*mod(k-1,ncol) 1-(stepy+hight)*(ceil(k/ncol)) width hight],...
            'Tag',sprintf('changelater%d',k),...
            'xtick',[],'ytick',[],'xlim',[mnx mxx],'ylim',[mny mxy],...
            'Parent',handles.hfeatures,'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
        %         axes(handles.hfaxes(k));
        h=title(sprintf('%s vs %s',handles.feature_names{i},handles.feature_names{j}));
        set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
            'Parent',handles.hfaxes(k));
        box on
        k=k+1;
    end
end

function closefeatures_callback(source,eventdata);
%does not really close, makes it invisible
handles=guidata(get(source,'UserData'));
set(handles.hfeatures,'Visible','Off');
set(handles.mainfig,'Visible','On');

function copy_to_new_window(source,event);
figure(100);
cla;
copyobj(get(source,'children'),gca);
%,'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add'