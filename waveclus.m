function waveclus
%this version does not use setappdata ot getappdtata, all data is
%transfered through handles
%create parameters, make sure that they do not overlap with existing ones

handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
handles.MINISI=[0.01:0.01:0.09 0.1:0.1:3 4:64 100:100:1000];
handles.MAX_CLUS=15;
handles.classify_space='features';
handles.classify_method= 'linear';
clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0 0 0.5;0 0 0];
set(0,'DefaultAxesColorOrder',clus_colors);
handles.colors= clus_colors;

%order is important
handles=create_mainfig(handles);
handles=create_supplfig(handles);
handles=create_supplfig2(handles);
handles.spikeaxes=[handles.spikeaxes handles.axesClust0];
%handles.isiaxes=[handles.isiaxes handles.axesISI0];

handles=wc_create_detailsfigure(handles);

figure(handles.mainfig);
% Update handles structure
guidata(handles.mainfig, handles);

function loadbutton_Callback(source, eventdata)
handles=guidata(get(source,'UserData'));
handles=wc_clean_handles(handles);

% Update handles structure
guidata(handles.mainfig, handles);
handles=guidata(get(source,'UserData'));
mainfolder=['Y:' filesep 'Data' filesep 'Sortcodes' filesep];
if exist(mainfolder,'dir')
    [filename, pathname] = uigetfile('*spikes_ch*.mat','Select file',mainfolder);
else
    [filename, pathname] = uigetfile('*spikes_ch*.mat','Select file');
end
if isequal([filename,pathname],[0,0]), return; end
%cd(pathname);
switch get(handles.htoread,'Value'),
    case 1, %ION
        t1=strfind(filename,'_ch');%minus is used to divide bname from channel
        t2=strfind(filename,'.mat');
        handles.bname = filename(1:t1-1);
        handles.thrsign = filename(t1+7:t2-1);
        handles.channel=str2double(filename(t1+3:t1+5));
        handles.filename=(filename(1:t2(1)-1));
    case 2, %UCLA
end
set(handles.textStatus,'string',sprintf('Loading %s',handles.filename),'fontsize',14);

q=load(sprintf('%s.mat',[pathname filesep handles.filename]));

handles.WC=q.par;
handles.WC.thr = q.thr;
handles.index=q.cluster_class(:,2);
handles.nspk=length(handles.index);
handles.ncl=max(q.cluster_class(:,1));
handles.nfeatures=size(q.features,2);


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


handles.sp_time=-handles.WC.w_pre*handles.WC.int_factor+1:1:handles.WC.w_post*handles.WC.int_factor;
handles.classtemp=q.classtemp;
handles.tree=q.tree;
handles.features=q.features;
handles.feature_names=q.feature_names;

%assign clusters as they were saved
for i=1:handles.ncl, handles.classind{i}=find(q.cluster_class(:,1)==i)'; end
%cluster zero
handles.classind{end+1}=setdiff(1:handles.nspk,[handles.classind{:}]);
clear q;

load([pathname filesep 'concatenation_info.mat'],'blocksamplesperchannel','wheretofindwhat','whattofindwhere','channels_to_process','sr');
us_idx=strfind(filename,'_');
n_file=str2double(filename(us_idx(2)+1:us_idx(3)-1));
block=whattofindwhere{handles.channel}{n_file}(1);
us_idx=strfind(pathname,filesep);
tens_fname=[pathname(1:us_idx(end-1)) 'WC_Block-' num2str(block) filesep 'datafilt_ch' sprintf('%03d.mat',handles.channel)];
if ~exist(tens_fname,'file'),
    handles.ts=[0 0];
    handles.ts_time=[0 1];
else
    q=load(tens_fname);
    handles.ts=double(q.data(1:round(10*handles.WC.sr)))*handles.WC.transform_factor;
    handles.ts_time=(1:length(handles.ts))/handles.WC.sr;
    clear q;
    set(handles.textStatus,'string',sprintf('Plotting %s',handles.filename));
    wc_plot_raw(handles);
end

handles=wc_plot_temperature(handles);
handles=wc_plot_spikes_and_ISI(handles);

set(handles.textStatus,'string',sprintf('%s',handles.filename));

% Update handles structure
guidata(handles.mainfig, handles);


function detailsbuttons_Callback(source,~)
% Plot details of selected cluster, a new figure
tag=get(source,'Tag');
i=str2double(tag(8:end));
handles=guidata(get(source,'UserData'));
handles.details_currcluster=i;

handles=wc_plot_details(handles);

% Update handles structure
guidata(handles.mainfig, handles);

function classifybutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));

%think about automatically plotting a new (updated) version of features
% %make feature plot invisible
% if isfield(handles,'hfeatures'), set(handles.hfeatures,'Visible','Off'); end

if get(handles.hclassify,'value') == 1,
    handles=wc_classifyrest(handles,2);
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

handles=wc_plot_spikes_and_ISI(handles);
if ~all(handles.ts==0), wc_plot_raw(handles); end
% Update handles structure
guidata(handles.mainfig, handles);


function tempmatchbutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));

%think about automatically plotting a new (updated) version of features
% %make feature plot invisible
% if isfield(handles,'hfeatures'), set(handles.hfeatures,'Visible','Off'); end

if get(handles.hclassify,'value') == 1,
    handles=wc_classifyrest(handles,1);
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

handles=wc_plot_spikes_and_ISI(handles);
if ~all(handles.ts==0), wc_plot_raw(handles); end
% Update handles structure
guidata(handles.mainfig, handles);

function changetempbutton_Callback(source, ~)
[temp min_clus]= ginput(1);                  %gets the mouse input
handles=guidata(get(source,'UserData'));
temp = round(temp)+1;
if temp < 1; temp=1;end                 %temp should be within the limits
if temp > handles.WC.num_temp; temp=handles.WC.num_temp; end
min_clus = round(min_clus);
handles.hhor=handles.hhor(1);
set(handles.hhor,'ydata',[min_clus min_clus]);
set(handles.hver,'xdata',[temp-1 temp-1]);

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

set(handles.hfix(fixed),'Value',0);
handles.fixed(fixed)=0;

tree=handles.tree;
toplot=find(tree(temp,5:end)>=min_clus);%find clusters which are big enough to plot


fixed_temps=handles.WC.clus_per_temp(:,fixed);
for c=1:numel(handles.hcol)
    delete(handles.hcol(c));
end


nnewclust=0;
new_temps=[];
for i=1:length(toplot),
    t=handles.classtemp{temp,toplot(i)};
    %remove spikes which belong to fixed clusters
    t=setdiff(t,fixed_classind_all);
    if ~isempty(t),
        nnewclust=nnewclust+1;
        handles.classind{nnewclust}=t;
        new_temps=[new_temps [temp;toplot(i)]];
    end
end

handles.WC.clus_per_temp=[new_temps fixed_temps];

handles.hcol=[];
colidx=1;
for c=1:size(handles.WC.clus_per_temp,2)
    x=handles.WC.clus_per_temp(1,c)-1;
    handles.hver(c)=line([x x],[1 max(max(tree(:,5:end)))*1.1],'linestyle',':','color','k');
    cc=handles.WC.clus_per_temp(2,c);
    y=tree(x+1,cc+4);
    handles.hcol(colidx)=scatter(x,y,50,handles.colors(colidx,:),'o','filled');
    colidx=colidx+1;
end

handles.ncl=nfixed+nnewclust;

% put fixed clusters to the end of all clusters
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

%cluster zero
handles.classind{end+1}=setdiff(1:handles.nspk,[handles.classind{:}]);
handles=wc_plot_spikes_and_ISI(handles);
if ~all(handles.ts==0), wc_plot_raw(handles); end
guidata(handles.mainfig, handles);

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
handles.classind=handles.classind(setdiff(1:length(handles.classind),i));

handles.ncl=handles.ncl-1;
handles.rejected=i;
handles=wc_plot_spikes_and_ISI(handles);
handles.rejected=1;
if ~all(handles.ts==0), wc_plot_raw(handles); end
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
handles=wc_plot_features_vs_features(handles);

function plottimecoursebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
handles=wc_plot_features_vs_time(handles);

function toplot_Callback(source,~)
handles=guidata(get(source,'UserData'));
%if data already loaded
if isfield(handles,'spikes'), handles=wc_plot_spikes_and_ISI(handles); end
% Update handles structure
guidata(handles.mainfig, handles);

function toread_Callback(source,eventdata)
% handles=guidata(get(source,'UserData'));
loadbutton_Callback(source, eventdata)

function savebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Saving %s',handles.filename));
handles=wc_saveresults(handles);
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
handles.isaGUI=1;
handles.mainfig=figure('Visible','on','Units','Normalized','Position',[0 0 1,0.9],...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Name','Wave_clus','NumberTitle','off',...
    'Tag','Mainfig');

figure(handles.mainfig);

ncol=5;
nrow=4;
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
handles.t=1;
handles.axesAllClusters=axes('position',[stepx 1-(stepy+hight)*( 2 ) width hight],...load
    'Tag','AllClusters','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.axesClust0=axes('position',[stepx+(width+stepx)*0 1-(stepy+hight)*( nrow ) width hight],...
    'Tag','Clust0','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.hclassify=uicontrol('units','normalized','Style','togglebutton','String','Classify','FontSize',12,...
    'Tag','classify',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.01 1-(stepy+hight)*( nrow )-stepy*0.9 0.05 0.035],...
    'Callback',{@classifybutton_Callback});
handles.htempmatch=uicontrol('units','normalized','Style','togglebutton','String','TempMatch','FontSize',12,...
    'Tag','Tempmatch',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.065 1-(stepy+hight)*( nrow )-stepy*0.9 0.05 0.035],...
    'Callback',{@tempmatchbutton_Callback});
handles.hexit=uicontrol('units','normalized','Style','pushbutton','String','Exit','FontSize',12,...
    'Tag','exit',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.12 1-(stepy+hight)*( nrow )-stepy*0.9 0.05 0.035],...
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
    'Position',[stepx 1-(stepy+hight)*( 3 )-stepy*0.65 0.15 0.03],...
    'Callback',{@changetempbutton_Callback});

i=0;
for r=2:nrow
    for c=1:ncol-1,
        i=i+1;
        handles.spikeaxes(i)=axes('position',[stepx+(width+stepx)*c 1-(stepy+hight)*r width hight],...
            'Tag',sprintf('Clust%d',i),'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
        handles.hfix(i) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
            'Tag',sprintf('Fix%d',i),...
            'UserData',handles.mainfig,...
            'Position',[stepx+(width+stepx)*c 1-(stepy+hight)*r-stepy*0.45 0.035 0.02],...
            'Callback',{@fixbuttons_Callback});
        handles.hreject(i) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
            'Tag',sprintf('Reject%d',i),...
            'UserData',handles.mainfig,...
            'Position',[stepx+(width+stepx)*c+0.04 1-(stepy+hight)*r-stepy*0.45 0.02 0.02],...
            'Callback',{@rejectbuttons_Callback});
        handles.hdetails(i) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
            'Tag',sprintf('Details%d',i),...
            'UserData',handles.mainfig,...
            'Position',[stepx+(width+stepx)*c+0.1 1-(stepy+hight)*r-stepy*0.45 0.05 0.02],...
            'Callback',{@detailsbuttons_Callback});
        handles.hclustergroup{i}=[handles.spikeaxes(i) handles.hfix(i) handles.hreject(i) handles.hdetails(i)];
        set(handles.hclustergroup{i},'Visible','Off');
    end
end

function handles=create_supplfig(handles)

handles.hsuppl = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig');
nrow=4;ncol=5;
stepx=0.04; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;
for i=1:ncol*2,
    j=i+12;
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
    handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j)];
    set(handles.hclustergroup{j},'Visible','Off');
end

function handles=create_supplfig2(handles)
handles.hsuppl2 = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters 2','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig2');
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
    handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j)];
    set(handles.hclustergroup{j},'Visible','Off');
end

function copy_to_new_window(source,~)
figure(100);
cla;
copyobj(get(source,'children'),gca);
