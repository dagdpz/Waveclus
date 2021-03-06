function handles=create_timecoursefig(handles);
if isfield(handles,'htimecourse'), 
    figure(handles.htimecourse);clf;
else
    if ishandle(handles.mainfig) %no gui opend
        handles.htimecourse=figure('Visible','Off','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
        'Name','Features vs time','NumberTitle','off','Color','w',...
        'UserData',handles.mainfig,...
        'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
        'CloseRequestFcn',{@closefeatures_callback},...
        'Tag','Features');
    else
         handles.htimecourse=figure('Visible','Off','Units','Normalized','Position',[0.01,0.01,0.9,0.9],...
        'Name','Features vs time','NumberTitle','off','Color','w',...
        'UserData',handles.mainfig,...
        'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1600 1080],...
        'Tag','Features');
    end
   
%     'Paperunits','Normalized','Paperorientation','Landscape','PaperPosition',[0.01 0.01 0.98 0.98],...
       
end

nf=handles.nfeatures;
nrow=ceil(sqrt(nf));
ncol=ceil(sqrt(nf));

stepx=0.01; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.01; hight=(1-(nrow+1)*stepy)/nrow;
        mny=min(min(handles.index));
        mxy=max(max(handles.index));
for i=1:nf,
    %for j=i+1:nf,
        mnx=min(min(handles.features(:,i)));
        mxx=max(max(handles.features(:,i)));
        handles.htaxes(i)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*(ceil(i/ncol)) width hight],...
            'Tag',sprintf('changelater%d',i),...
            'xtick',[],'ytick',[],'xlim',[mnx mxx],'ylim',[mny mxy],...
            'Parent',handles.htimecourse,'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
%         axes(handles.hfaxes(k));
        h=title(sprintf('time vs %s',handles.feature_names{i}));
        set(h,'Units','Normalized','Position',[0.01 .01 0],'verticalalignment','bottom','horizontalalignment','left',...
            'Parent',handles.htaxes(i));
        box on
    %end    
end


function closefeatures_callback(source,eventdata);
%does not really close, makes it invisible
handles=guidata(get(source,'UserData'));
set(handles.htimecourse,'Visible','Off');
set(handles.mainfig,'Visible','On');

function copy_to_new_window(source,event);
figure(100);
cla;
copyobj(get(source,'children'),gca);
%,'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add'