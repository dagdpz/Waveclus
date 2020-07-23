function handles=wc_create_featuresfig(handles)

% if isfield(handles,'hfeatures') && isobject(handles.hfeatures),
%     figure(handles.hfeatures);clf;
% else
    if ishandle(handles.mainfig)
        handles.hfeatures=figure('Visible','On','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
            'Name','Features','NumberTitle','off','Color','w',...
            'UserData',handles.mainfig,...
            'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...%            'CloseRequestFcn',{@closefeatures_callback},...
            'Tag','Features');
    else %no gui opend
        handles.hfeatures=figure('Visible','On','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
            'Name','Features','NumberTitle','off','Color','w',...
            'UserData',handles.mainfig,...
            'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
            'Tag','Features');
    end
%end

nf=handles.nfeatures;
np=nf;%*(nf-1)/2;
nrow=floor(sqrt(np));
ncol=ceil(sqrt(np));

stepx=0.01; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.01; hight=(1-(nrow+1)*stepy)/nrow;
k=2;
for i=2 %:nf,
    for j=i+1:nf,
    mnx=prctile(handles.features(:,i),0.1);
    mxx=prctile(handles.features(:,i),99.9);
    mny=prctile(handles.features(:,j),0.1);
    mxy=prctile(handles.features(:,j),99.9);
%         mnx=min(min(handles.features(:,i)));
%         mny=min(min(handles.features(:,j)));
%         mxx=max(max(handles.features(:,i)));
%         mxy=max(max(handles.features(:,j)));
        handles.hfaxes(k)=axes('position',[stepx+(width+stepx)*mod(k-1,ncol) 1-(stepy+hight)*(ceil(k/ncol)) width hight],...
            'Tag',sprintf('changelater%d',k),...
            'xtick',[],'ytick',[],'xlim',[mnx mxx]*1.1,'ylim',[mny mxy]*1.1,...
            'Parent',handles.hfeatures,'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
        h=title(sprintf('%s vs %s',handles.feature_names{j},handles.feature_names{i}));
        set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
            'Parent',handles.hfaxes(k));
        box on
        k=k+1;
    end
end

% function closefeatures_callback(source,eventdata)
% %does not really close, makes it invisible
% handles=guidata(get(source,'UserData'));
% set(handles.hfeatures,'Visible','Off');
% set(handles.mainfig,'Visible','On');

function copy_to_new_window(source,event)
figure(100);
cla;
copyobj(get(source,'children'),gca);