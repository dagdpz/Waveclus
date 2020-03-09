function wave_clus3new3

% set(groot,'DefaultFigureGraphicsSmoothing','off')
%this version does not use setappdata ot getappdtata, all data is
%transfered through handles
%create parameters, make sure that they do not overlap with existing ones

handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
handles.MINISI=[0.01:0.01:0.09 0.1:0.1:3 4:64 100:100:1000];
handles.MAX_CLUS=23;
% handles.classify_space='features';
% handles.classify_space='spikeshapes';
handles.classify_space='spikeshapesfeatures';
% handles.classify_method= 'diaglinear'; 
handles.classify_method= 'linear';

clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5;0 0 0];
set(0,'DefaultAxesColorOrder',clus_colors);
% handles.colors='brgcmybrgcmybrgcmybrgcmy';
handles.colors= clus_colors;

%order is important
handles=create_mainfig(handles);
handles=create_supplfig(handles);
handles=create_supplfig2(handles);
handles.spikeaxes=[handles.spikeaxes handles.axesClust0];
handles.isiaxes=[handles.isiaxes handles.axesISI0];

handles=create_detailsfigure(handles);

figure(handles.mainfig);
% Update handles structure
guidata(handles.mainfig, handles);

function loadbutton_Callback(source, eventdata)
handles=guidata(get(source,'UserData'));

handles=clean4new(handles);

% Update handles structure
guidata(handles.mainfig, handles);

% t=cputime;  fprintf('Start loading\n');

handles=guidata(get(source,'UserData'));
% t=cputime-t; fprintf('Got handles \t\t%1.0f\n',t); t=cputime;

[filename, pathname] = uigetfile('*spikes_ch*.mat','Select file');
if isequal([filename,pathname],[0,0]), return; end
cd(pathname);
switch get(handles.htoread,'Value'),
    case 1, %ION
        t1=strfind(filename,'_ch');%minus is used to divide bname from channel
        t2=strfind(filename,'.mat');
        handles.bname = filename(1:t1-1);
        handles.thrsign = filename(t1+7:t2-1);
        handles.channel=str2double(filename(t1+3:t1+5));
        handles.filename=(filename(1:t2(1)-1));
        
%         bname=filename(7:t1(1));%skip word times
%         handles.bname=bname;
%         ch=str2double(filename(t1(end)+1:t2(1)-1));
%         handles.channel=ch;
%         handles.filename=(filename(7:t2(1)-1));
    case 2, %UCLA
end
set(handles.textStatus,'string',sprintf('Loading %s',handles.filename),'fontsize',14);

q=load(sprintf('%s.mat',handles.filename));

%%temp
% q.features=11;
% q.par.min_clus

handles.par=q.par; 
handles.par.thr = q.thr;
handles.index=q.cluster_class(:,2);
handles.nspk=length(handles.index);
handles.ncl=max(q.cluster_class(:,1));
handles.nfeatures=size(q.features,2);
% --------
% handles.feature_names = q.feature_names;
%  handles.spikecomp = q.spikecomp;
% --------

handles.min_clus=q.par.min_clus;
if isfield(q.par,'temp'),handles.temp=q.par.temp; else handles.temp=0; end

if isfield(q,'spikes')
    handles.spikes=q.spikes;
else
    spikesfile=sprintf('spikes_%s.mat',handles.filename);
    qq=load(spikesfile);
    handles.spikes=qq.spikes;
    clear qq;
end   

% handles.sp_time=(-handles.par.w_pre+1:1:handles.par.w_post)/handles.par.sr*1000;%time scale for spike shapes in ms
%-----------------------------
handles.sp_time=-handles.par.w_pre*handles.par.int_factor+1:1:handles.par.w_post*handles.par.int_factor;
%----------------------------------
handles.classtemp=q.classtemp;
handles.tree=q.tree;
handles.features=q.features;
handles.feature_names=q.feature_names;

% t=cputime-t; fprintf('end assigning handles \t\t%1.0f\n',t); t=cputime;

%assign clusters as they were saved
for i=1:handles.ncl, handles.classind{i}=find(q.cluster_class(:,1)==i)'; end
%cluster zero
handles.classind{end+1}=setdiff(1:handles.nspk,[handles.classind{:}]);
clear q;

% t=cputime-t; fprintf('end assiging classind \t\t%1.0f\n',t); t=cputime;
load([pathname filesep 'concatenation_info.mat'],'blocksamplesperchannel','wheretofindwhat','whattofindwhere','channels_to_process','sr');
us_idx=strfind(filename,'_');
n_file=str2double(filename(us_idx(2)+1:us_idx(3)-1));
block=whattofindwhere{handles.channel}{n_file}(1);
us_idx=strfind(pathname,filesep);
tens_fname=[pathname(1:us_idx(end-1)) 'WC_Block-' num2str(block) filesep 'datafilt_ch' sprintf('%03d.mat',handles.channel)];
if ~exist(tens_fname,'file'), 
    %     error('10 seconds of data file does not exist');
    handles.ts=[0 0];
    handles.ts_time=[0 1];
else 
    q=load(tens_fname); 
    handles.ts=double(q.data(1:round(10*handles.par.sr)))*handles.par.transform_factor;
    handles.ts_time=(1:length(handles.ts))/handles.par.sr; 
    %thesholds from spikes file?
    clear q; 
    set(handles.textStatus,'string',sprintf('Plotting %s',handles.filename));
    % t=cputime-t; fprintf('end loading ts \t\t%1.0f\n',t); t=cputime;
    
    %plotting
    plot_ts(handles);
end

% t=cputime-t; fprintf('end plotting ts \t\t%1.0f\n',t); t=cputime;
handles=plot_temperature(handles);
% t=cputime-t; fprintf('end plotting temperatures \t\t%1.0f\n',t); t=cputime;
handles=plot_spikes_isi(handles);
% t=cputime-t; fprintf('end plotting spikes fig \t\t%1.0f\n',t); t=cputime;

set(handles.textStatus,'string',sprintf('%s',handles.filename));

% %placed here since it does depend on loaded data (number of features)
% if isfield(handles,'hfeatures'), delete(handles.hfeatures); end
% handles=create_featuresfig(handles);

% t=cputime-t; fprintf('end creating features fig \t\t%1.0f\n',t); t=cputime;

% Update handles structure
guidata(handles.mainfig, handles);
% t=cputime-t; fprintf('End updating handles \t\t%1.0f\n',t); t=cputime;


function detailsbuttons_Callback(source,~)
% Plot details of selected cluster, a new figure
tag=get(source,'Tag');
i=str2double(tag(8:end));
handles=guidata(get(source,'UserData'));
handles.details_currcluster=i;

handles=plot_details(handles);

% Update handles structure
guidata(handles.mainfig, handles);

function classifybutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));

%think about automatically plotting a new (updated) version of features

% %make feature plot invisible
% if isfield(handles,'hfeatures'), set(handles.hfeatures,'Visible','Off'); end

if get(handles.hclassify,'value') == 1,
    handles=classifyrest3(handles);
    %remove fix from forced clusters
    for i=find(handles.forced),
        set(handles.hfix(i),'Value',0);
        handles.fixed(i)=0;
    end
    set(handles.hclassify,'String','Classified');
else
    for i=1:handles.ncl,
        if ~handles.fixed(i), handles.classind{i}=handles.classind_unforced{i};
        else handles.fixed(i)=0; set(handles.hfix(i),'Value',0); end
    end
    handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);    
    handles.forced(find(handles.forced))=0;
    set(handles.hclassify,'String','Classify');
end

handles=plot_spikes_isi(handles);
if ~all(handles.ts==0), plot_ts(handles); end
% Update handles structure
guidata(handles.mainfig, handles);


function tempmatchbutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));

if get(handles.htempmatch,'value') == 1,
    handles=tempmatch(handles);
    for i=find(handles.forced),
        set(handles.hfix(i),'Value',0);
        handles.fixed(i)=0;
    end
    set(handles.htempmatch,'String','Classified');
% else
%     for i=1:handles.ncl,
%         if ~handles.fixed(i), handles.classind{i}=handles.classind_unforced{i};
%         else handles.fixed(i)=0; set(handles.hfix(i),'Value',0); end
%     end
%     handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);    
%     handles.forced(find(handles.forced))=0;
%     set(handles.htempmatch,'String','Classify');
end
        
handles=plot_spikes_isi(handles);
if ~all(handles.ts==0), plot_ts(handles); end
% Update handles structure
guidata(handles.mainfig, handles);



function changetempbutton_Callback(source, ~)
% ti=cputime; fprintf('Start change temperature \t\t%1.0f\n',ti); ti=cputime;
[temp min_clus]= ginput(1);                  %gets the mouse input
% ti=cputime-ti; fprintf('End input \t\t%1.0f\n',ti); ti=cputime;

handles=guidata(get(source,'UserData'));

temp = round(temp);
if temp < 1; temp=1;end                 %temp should be within the limits
if temp > handles.par.num_temp; temp=handles.par.num_temp; end
min_clus = round(min_clus);
set(handles.hhor,'ydata',[min_clus min_clus]);
set(handles.hver,'xdata',[temp temp]);

% ti=cputime-ti; fprintf('Got handles \t\t%1.0f\n',ti); ti=cputime;

% axes(handles.axesTemp)

% %make feature plot invisible
% if isfield(handles,'hfeatures'), set(handles.hfeatures,'Visible','Off'); end
%reset forced
handles.forced(find(handles.forced))=0;
set(handles.hclassify,'String','Classify');
set(handles.hclassify,'Value',0);

%--------
set(handles.htempmatch,'String','TempMatch');
set(handles.htempmatch,'Value',0);
%-------
    
handles.classind={};
handles.temp=temp;
handles.min_clus=min_clus;

fixed=find(handles.fixed);%all fixed clusters
nfixed=length(fixed);
fixed_classind_all=[handles.fixed_classind{fixed}];%all spike indices in fixed cluster
    
toplot=find(handles.tree(temp,5:end)>=min_clus);%find clusters which are big enough to plot
nnewclust=0;
for i=1:length(toplot),
    t=handles.classtemp{temp,toplot(i)};
    %remove spikes which belong to fixed clusters
    t=setdiff(t,fixed_classind_all);
    if ~isempty(t), 
        nnewclust=nnewclust+1;
        handles.classind{nnewclust}=t;
    end
end
handles.ncl=nfixed+nnewclust;

set(handles.hfix(fixed),'Value',0);
handles.fixed(fixed)=0;
%put fixed clusters to the end of all clusters
for i=nnewclust+1:handles.ncl,
    j=i-nnewclust;
    handles.classind{i}=handles.fixed_classind{fixed(j)};
    handles.fixed_classind{fixed(j)}=[];
    set(handles.hfix(i),'Value',1);
    handles.fixed(i)=1;
end
for i=nnewclust+1:handles.ncl,
    handles.fixed_classind{i}=handles.classind{i};
end

% %put fixed clusters to the end of all unfixed clusters, start from the end
% %to avoid putting one and the same cluster into all fixed
% for i=k:handles.ncl,
%     j=i-k+1;
%     handles.classind{i}=handles.fixed_classind{fixed(j)};
%     handles.fixed_classind{fixed(j)}=[];
%     
%     handles.fixed_classind{ii}=handles.classind{ii};%move fixed clusters onto different plot
%     handles.fixed(fixed(j))=0;
%     handles.fixed(ii)=1;
%     set(handles.hfix(fixed(j)),'Value',0);
%     set(handles.hfix(ii),'Value',1);
% end
%cluster zero
handles.classind{end+1}=setdiff(1:handles.nspk,[handles.classind{:}]);
% ti=cputime-ti; fprintf('end assinging clusters \t\t%1.0f\n',ti); ti=cputime;
handles=plot_spikes_isi(handles);
if ~all(handles.ts==0), plot_ts(handles); end
% ti=cputime-ti; fprintf('end plotting spikes \t\t%1.0f\n',ti); ti=cputime;
% Update handles structure
guidata(handles.mainfig, handles);
% ti=cputime-ti; fprintf('End updating handles \t\t%1.0f\n',ti); ti=cputime;

function rejectbuttons_Callback(source,~)
% Reject cluster
tag=get(source,'Tag');
i=str2double(tag(7:end));

handles=guidata(get(source,'UserData'));
%reset forced
handles.forced(find(handles.forced))=0;
set(handles.hclassify,'String','Classify');
set(handles.hclassify,'Value',0);

if find(handles.fixed),
    handles.fixed(i:end-1)=handles.fixed(i+1:end);
    handles.fixed(end)=0;
    handles.fixed_classind={handles.fixed_classind{setdiff(1:length(handles.fixed_classind),i)}};
    handles.fixed_classind{end+1}=[];
    for j=i:length(handles.hfix)-1,
        set(handles.hfix(j),'Value',get(handles.hfix(j+1),'Value'));
    end
    set(handles.hfix(end),'Value',0);
end
    
handles.classind{end}=[handles.classind{end} handles.classind{i}];
handles.classind={handles.classind{setdiff(1:length(handles.classind),i)}};

handles.ncl=handles.ncl-1;
handles.rejected=i;
handles=plot_spikes_isi(handles);
handles.rejected=1;
if ~all(handles.ts==0), plot_ts(handles); end
% Update handles structure
guidata(handles.mainfig, handles);

function fixbuttons_Callback(source,~)
tag=get(source,'Tag');
i=str2double(tag(4:end));
handles=guidata(get(source,'UserData'));
switch get(source,'value');
    case 1, handles.fixed(i)=1;handles.fixed_classind{i}=handles.classind{i};
    case 0, handles.fixed(i)=0;handles.fixed_classind{i}=[];
end
% Update handles structure
guidata(handles.mainfig, handles);


function plotfeaturesbutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
handles=plot_features2(handles);

function plottimecoursebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
%handles=plot_features2(handles);
handles=plot_timecourse(handles);
% set(handles.mainfig,'Visible','off')

function toplot_Callback(source,~)
handles=guidata(get(source,'UserData'));
%if data already loaded
if isfield(handles,'spikes'), handles=plot_spikes_isi(handles); end
% Update handles structure
guidata(handles.mainfig, handles);

function toread_Callback(source,eventdata)
% handles=guidata(get(source,'UserData'));
loadbutton_Callback(source, eventdata)

function savebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Saving %s',handles.filename));
handles=saveresults2(handles);
figure(handles.mainfig);
set(handles.textStatus,'string',sprintf('Saved %s',handles.filename));

%append only changed variables, classind, handles (temp, min_clus),
%unforced version?

function exitbutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
eval('delete(handles.hfeatures)','');
eval('delete(handles.hsuppl)','');
eval('delete(handles.hsuppl2)','');
eval('delete(handles.hdetailsfig)','');
eval('delete(handles.mainfig)','');

function handles=create_mainfig(handles)
handles.mainfig=figure('Visible','on','Units','Normalized','Position',[0 0 1,0.9],...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Name','Wave_clus','NumberTitle','off',...
    'Tag','Mainfig');
% 'Paperunits','Normalized','Paperorientation','portrait','PaperPosition',[0.01 0.01 0.98 0.98],...

figure(handles.mainfig);
nrow=3;ncol=5;
stepx=0.02; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;

handles.htoread = uicontrol('units','normalized','Style','popupmenu','FontSize',12,...
    'String',{'ION','UCLA'},...
    'UserData',handles.mainfig,...
    'Position',[stepx 1-stepy*0.9 0.15 0.035],...
    'Callback',{@toread_Callback});
handles.hload=uicontrol('units','normalized','Style','pushbutton','String','Load','FontSize',12,...
    'Tag','load',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*1 1-stepy*0.9 0.05 0.035],...
    'Callback',{@loadbutton_Callback});
handles.textStatus=uicontrol('units','normalized','Style','Text','String','','FontSize',12,...
    'Tag','textstatus',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*2 1-stepy*0.9 0.35 0.035]);
handles.hsave=uicontrol('units','normalized','Style','pushbutton','String','Save','FontSize',12,...
    'Tag','save',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*4 1-stepy*0.9 0.05 0.035],...
    'Callback',{@savebutton_Callback});

handles.axesTS=axes('position',[stepx 1-(stepy+hight)*( 1 )+stepy*0.5 1-2*stepx hight*0.9],...
    'Tag','TS',...
    'UserData',handles.mainfig,...
    'ButtonDownFcn',{@copy_to_new_window});
% handles.t=1;;
handles.axesAllClusters=axes('position',[stepx 1-(stepy+hight)*( 2 ) width hight],...
    'Tag','AllClusters','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.axesClust0=axes('position',[stepx+(width+stepx)*4 1-(stepy+hight)*( 2 ) width hight],...
    'Tag','Clust0','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.axesISI0=axes('position',[stepx+(width+stepx)*4 1-(stepy+hight)*( 3 )+hight*0.1 width hight*0.9],...
    'Tag','ISI0','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.hclassify=uicontrol('units','normalized','Style','togglebutton','String','Classify','FontSize',12,...
    'Tag','classify',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*4+0.01 1-(stepy+hight)*( 3 )-stepy*0.9 0.05 0.035],...
    'Callback',{@classifybutton_Callback});
handles.htempmatch=uicontrol('units','normalized','Style','togglebutton','String','TempMatch','FontSize',12,...
    'Tag','Tempmatch',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*4+0.065 1-(stepy+hight)*( 3 )-stepy*0.9 0.05 0.035],...
    'Callback',{@tempmatchbutton_Callback});
handles.hexit=uicontrol('units','normalized','Style','pushbutton','String','Exit','FontSize',12,...
    'Tag','exit',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*4+0.12 1-(stepy+hight)*( 3 )-stepy*0.9 0.05 0.035],...
    'Callback',{@exitbutton_Callback});

handles.axesTemp=axes('position',[stepx 1-(stepy+hight)*( 3 )+0.1*hight width hight],...
        'Tag','Temperature');
handles.hplotfeatures=uicontrol('units','normalized','Style','pushbutton','String','Features','FontSize',12,...
    'Tag','plotfeatures',...
    'UserData',handles.mainfig,...
    'Position',[stepx 1-(stepy+hight)*( 1 )-stepy*0.9 0.05 0.03],...
    'Callback',{@plotfeaturesbutton_Callback});
handles.hplotfeatures=uicontrol('units','normalized','Style','pushbutton','String','Timecourse','FontSize',12,...
    'Tag','plottimecourse',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.05 1-(stepy+hight)*( 1 )-stepy*0.9 0.05 0.03],...
    'Callback',{@plottimecoursebutton_Callback});
handles.htoplot = uicontrol('units','normalized','Style','popupmenu','FontSize',12,...
    'String',{'All','Average','5%','10%','20%','50%'},...
    'Value',3,...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.12 1-(stepy+hight)*( 1 )-stepy*0.9 0.06 0.03],...
    'Callback',{@toplot_Callback});

handles.hchangetemp=uicontrol('units','normalized','Style','pushbutton','String','Change Temp','FontSize',12,...
    'Tag','ChangeTemp',...
    'UserData',handles.mainfig,...
    'Position',[stepx 1-(stepy+hight)*( 3 )-stepy*0.75 0.15 0.03],...
    'Callback',{@changetempbutton_Callback});

for i=1:3,
    j=i+1+ncol;

    handles.spikeaxes(i)=axes('position',[stepx+(width+stepx)*mod(j-1,ncol) 1-(stepy+hight)*( 2 ) width hight],...
        'Tag',sprintf('Clust%d',i),'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');

    handles.hfix(i) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
        'Tag',sprintf('Fix%d',i),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(j-1,ncol) 1-(stepy+hight)*( 3 )-stepy*0.75 0.035 0.02],...
        'Callback',{@fixbuttons_Callback});
    handles.hreject(i) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
        'Tag',sprintf('Reject%d',i),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(j-1,ncol)+0.04 1-(stepy+hight)*( 3 )-stepy*0.75 0.02 0.02],...
        'Callback',{@rejectbuttons_Callback});
    handles.hdetails(i) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
        'Tag',sprintf('Details%d',i),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(j-1,ncol)+0.1 1-(stepy+hight)*( 3 )-stepy*0.75 0.05 0.02],...
        'Callback',{@detailsbuttons_Callback});
    
    handles.isiaxes(i)=axes('position',[stepx+(width+stepx)*mod(j-1,ncol) 1-(stepy+hight)*( 3 ) width hight],...
        'Tag',sprintf('ISI%d',i),'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
    
    handles.hclustergroup{i}=[handles.spikeaxes(i) handles.hfix(i) handles.hreject(i) handles.hdetails(i) handles.isiaxes(i)];
    set(handles.hclustergroup{i},'Visible','Off');
end
% TO DO, add context menu to define the type of unit
% handles.hcmenu=uicontextmenu;
% labels={'PTN','SU','MU','nSU','noise'};
% for i=1:length(labels),
%     handles.hm(i)=uimenu(handles.hcmenu,'Label',labels{i});
% end
% set(handles.mainfig,'UIContextMenu',handles.hcmenu);

function handles=create_supplfig(handles)

handles.hsuppl = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig');
% 'Paperunits','Normalized','Paperorientation','portrait','PaperPosition',[0.01 0.01 0.98 0.98],...
nrow=4;ncol=5;
stepx=0.04; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;
for i=1:ncol*2,
    j=i+3;%first three on the first plot
    handles.spikeaxes(j)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+1 ) width hight],...
        'Tag',sprintf('Clust%d',j));

    handles.hfix(j) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
        'Tag',sprintf('Fix%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.035 0.02],...
        'Callback',{@fixbuttons_Callback});
    handles.hreject(j) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
        'Tag',sprintf('Reject%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol)+0.04 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.02 0.02],...
        'Callback',{@rejectbuttons_Callback});
    handles.hdetails(j) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
        'Tag',sprintf('Details%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol)+0.1 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.05 0.02],...
        'Callback',{@detailsbuttons_Callback});
    
    handles.isiaxes(j)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+2 ) width hight],...
        'Tag',sprintf('ISI%d',j));
    handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j) handles.isiaxes(j)];
    set(handles.hclustergroup{j},'Visible','Off');
end

function handles=create_supplfig2(handles)
handles.hsuppl2 = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters 2','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig2');
% 'Paperunits','Normalized','Paperorientation','portrait','PaperPosition',[0.01 0.01 0.98 0.98],...
nrow=4;ncol=5;
stepx=0.04; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;
for i=1:ncol*2,
    j=i+3+ncol*2;%first three on the first plot
    handles.spikeaxes(j)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+1 ) width hight],...
        'Tag',sprintf('Clust%d',j));

    handles.hfix(j) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
        'Tag',sprintf('Fix%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.035 0.02],...
        'Callback',{@fixbuttons_Callback});
    handles.hreject(j) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
        'Tag',sprintf('Reject%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol)+0.04 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.02 0.02],...
        'Callback',{@rejectbuttons_Callback});
    handles.hdetails(j) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
        'Tag',sprintf('Details%d',j),...
        'UserData',handles.mainfig,...
        'Position',[stepx+(width+stepx)*mod(i-1,ncol)+0.1 1-(stepy+hight)*( (i>ncol)*2+2 )-stepy*0.75 0.05 0.02],...
        'Callback',{@detailsbuttons_Callback});
    
    handles.isiaxes(j)=axes('position',[stepx+(width+stepx)*mod(i-1,ncol) 1-(stepy+hight)*( (i>ncol)*2+2 ) width hight],...
        'Tag',sprintf('ISI%d',j));
    handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j) handles.isiaxes(j)];
    set(handles.hclustergroup{j},'Visible','Off');
end

% function testbuttondown(source,eventdata);
% handles=guidata(get(source,'UserData'));
% handles.t=~handles.t;
% fprintf('%d\n',handles.t);
% % Update handles structure
% guidata(handles.mainfig, handles);



function copy_to_new_window(source,~)
figure(100);
cla;
copyobj(get(source,'children'),gca);
% if isempty(h), return; end
% for i=1:length(h),
%     get()
% if ~isfield(h,'xdata'), return; end
% set('xdata',get(source,'xdata'),'ydata',get(source,'ydata'))
