function waveclus
%this version does not use setappdata ot getappdtata, all data is transfered through handles
%create parameters, make sure that they do not overlap with existing ones

handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
handles.MINISI=[0.01:0.01:0.09 0.1:0.1:3 4:64 100:100:1000];
handles.MAX_CLUS=15;
handles.classify_space='features';
handles.classify_method= 'linear';

handles.datafolder=['Y:' filesep 'Data' filesep 'Sortcodes' filesep];
clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.6 0 0; 0 0.6 0; 0 0 0.6; 0 0 0.6;0 0 0];
set(0,'DefaultAxesColorOrder',clus_colors);
handles.colors= clus_colors;

%order is important
handles=create_mainfig(handles);
handles=create_supplfig(handles);
handles=create_supplfig2(handles);
handles.spikeaxes=[handles.spikeaxes handles.axesClust0];
handles=wc_create_detailsfigure(handles);

figure(handles.mainfig);
% Update handles structure
guidata(handles.mainfig, handles);

%% file management buttons

function loadbutton_Callback(source, eventdata)
handles=guidata(get(source,'UserData'));
ts_time=handles.ts_time;
handles=wc_clean_handles(handles);
handles.ts_time=ts_time;
% Update handles structure
guidata(handles.mainfig, handles);
handles=guidata(get(source,'UserData'));

if isfield(handles,'pathname')
   currentfolder=handles.pathname;
else
   currentfolder=handles.datafolder;
end
if exist(currentfolder,'dir')
    [filename, pathname] = uigetfile('*spikes_ch*.mat','Select file',currentfolder);
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
        handles.pathname=pathname;
    case 2, %UCLA
end
set(handles.textStatus,'string',sprintf('Loading %s',handles.filename),'fontsize',14);

q=load(sprintf('%s.mat',[pathname filesep handles.filename]));
handles=handles_from_data(handles,q);

clear q;

% Update handles structure
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function handles=handles_from_data(handles,q)
handles.WC=q.par;
handles.WC.thr = q.thr;
handles.index=q.cluster_class(:,2);
handles.nspk=length(handles.index);
handles.ncl=max(q.cluster_class(:,1));
handles.nfeatures=size(q.features,2);
handles.plotted(1:handles.ncl)=1;

handles.min_clus=q.par.min_clus;
if isfield(q.par,'temp')
    handles.temp=q.par.temp; 
else
    handles.temp=0;
end
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


handles=wc_plot_temperature(handles);
guidata(handles.mainfig, handles);
plot_timecourse(handles.plot_timecourse);

function nextbutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
timesfile=sprintf('%s.mat',handles.filename);
allfiles=dir([handles.pathname filesep 'dataspikes*.mat']);
allfiles={allfiles.name};
idx=find(strcmp(allfiles,timesfile));

fileinvalid=1;
while fileinvalid
idx=mod(idx,numel(allfiles))+1;
handles.filename=allfiles{idx}(1:end-4);

ts_time=handles.ts_time;
handles=wc_clean_handles(handles);
handles.ts_time=ts_time;

set(handles.textStatus,'string',sprintf('Loading %s',handles.filename),'fontsize',14);

q=load(sprintf('%s.mat',[handles.pathname filesep handles.filename]));
if isfield(q,'features')
    fileinvalid=0;
end
end
if isfield(handles,'spike_indicators')
handles=rmfield(handles,'spike_indicators');
end

handles=handles_from_data(handles,q);
clear q;
guidata(handles.mainfig, handles); % Update handles structure!

function savebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Saving %s',handles.filename));
handles=guidata(get(source,'UserData'));
handles=wc_plot_features_vs_features(handles);
handles=guidata(get(source,'UserData'));
handles=wc_plot_features_vs_time(handles);
handles=wc_saveresults(handles);
close(handles.htimecourse);
close(handles.hfeatures)
figure(handles.mainfig);
set(handles.textStatus,'string',sprintf('Saved %s',handles.filename));

function toread_Callback(source,eventdata)
loadbutton_Callback(source, eventdata)

function detailsbuttons_Callback(source,~) %not used?
% Plot details of selected cluster, a new figure
tag=get(source,'Tag');
i=str2double(tag(8:end));
handles=guidata(get(source,'UserData'));
handles.details_currcluster=i;
handles=wc_plot_details(handles);
guidata(handles.mainfig, handles); % Update handles structure

%% Temperature & Classify

function classifybutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
do_or_undo=get(source,'Value');
accept_classification(handles);
if do_or_undo == 1,
    method=get_classify_version_from_tag(source);
    handles=wc_classifyrest(handles,method);
    set(source,'String','Undo');
    set(source,'Value',1);
else
    for i=1:handles.ncl,
        if ~handles.fixed(i), handles.classind{i}=intersect(handles.classind{i},[handles.classind_unforced{1:end-1}]); end
    end
    handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);
end
handles=wc_plot_spikes_and_ISI(handles);
if ~all(handles.ts==0), wc_plot_raw(handles); end
guidata(handles.mainfig, handles); % Update handles structure
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function accept_classification(handles)
set(handles.hclassify,'String','Class T');
set(handles.hclassify,'Value',0);
set(handles.htempmatch,'String','Classify');
set(handles.htempmatch,'Value',0);
set(handles.htempmatch2,'String','Near');
set(handles.htempmatch2,'Value',0);
set(handles.htempmatch3,'String','NearT');
set(handles.htempmatch3,'Value',0);
guidata(handles.mainfig, handles);

function method=get_classify_version_from_tag(source)
tag=get(source,'Tag');
switch tag
    case 'Class T'
        method=1;
    case 'Classify'
        method=2;
    case 'Near'
        method=3;
    case 'NearT'
        method=4;
end

function changetempbutton_Callback(source, ~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
%% Get the mouse input
[temp, min_clus, buttn]= ginput(1);                  
if buttn==3
    return;
end
temp = round(temp)+1;
if temp < 1; temp=1;end                 %temp should be discrete and within the limits
if temp > handles.WC.num_temp; temp=handles.WC.num_temp; end
min_clus = round(min_clus);
handles.temp        =temp;
handles.min_clus    =min_clus;
handles.hhor        =handles.hhor(1);
handles.classind    ={};
set(handles.hhor,'ydata',[min_clus min_clus]);
set(handles.hver,'xdata',[temp-1 temp-1]);

%% Treat fixed : Limit to only valid clusters
handles.fixed   =handles.fixed(1:handles.ncl);
handles.hcol    =handles.hcol(1:handles.ncl);
fixed=find(handles.fixed);%all fixed clusters
nfixed=length(fixed);
fixed_classind_all=[handles.fixed_classind{fixed}];%all spike indices in fixed cluster

%deleting all unit markers from temperature plot
for c=1:numel(handles.hcol)
    if ishghandle(handles.hcol(c)) 
    delete(handles.hcol(c));
    end
end
handles.hcol=[];

%% get new clusters
tree=handles.tree;
toplot=find(tree(temp,5:end-1)>=min_clus);%find clusters which are big enough to plot (-1 temporary debug)

nnewclust=0;
new_temps=[];
for i=1:length(toplot),
    t=handles.classtemp{temp,toplot(i)};    
    t=setdiff(t,fixed_classind_all); %remove spikes which belong to fixed clusters
    if ~isempty(t),
        nnewclust=nnewclust+1;
        handles.classind{nnewclust}=t;
        new_temps=[new_temps [temp;toplot(i)]];
    end
end

% cluster markers
fixed_temps=handles.WC.clus_per_temp(:,fixed(fixed<=size(handles.WC.clus_per_temp,2))); %fixed shouldnt be bigger than number of units with marker though
handles.WC.clus_per_temp=[new_temps fixed_temps];

%% draw vertical line and cluster markers
handles.hcol=[];
colidx=1;
for c=1:size(handles.WC.clus_per_temp,2)
    x=handles.WC.clus_per_temp(1,c)-1;
    handles.hver(c)=line([x x],[1 max(max(tree(:,5:end)))*1.1],'linestyle',':','color','k');
    cc=handles.WC.clus_per_temp(2,c);
    y=tree(x+1,cc+4);
    handles.hcol(colidx)=scatter(x,y,50,handles.colors(mod(colidx-1,size(handles.colors,1))+1,:),'o','filled');
    colidx=colidx+1;
end

% put fixed clusters to the end of all clusters and arrange the fix checkbox
handles.ncl=nfixed+nnewclust;
for i=1:nnewclust
    set(handles.hfix(i),'Value',0);
    handles.fixed(i)=0;
end
for i=nnewclust+1:handles.ncl,
    j=i-nnewclust;
    handles.classind{i}=handles.fixed_classind{fixed(j)};
    handles.fixed_classind{fixed(j)}=[];
    set(handles.hfix(i),'Value',1);
    handles.fixed(i)=1;
end
%cluster zero
handles.classind{end+1}=setdiff(1:handles.nspk,[handles.classind{:}]);

% dont think we need this
for i=nnewclust+1:handles.ncl,
    handles.fixed_classind{i}=handles.classind{i};
end

handles=wc_plot_spikes_and_ISI(handles);
if ~all(handles.ts==0), wc_plot_raw(handles); end
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

%% Feature plots 

function plotfeaturesbutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
wc_plot_features_vs_features(handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function plottimecoursebutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
wc_plot_features_vs_time(handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function toplot_Callback(source,~)
handles=guidata(get(source,'UserData'));
%if data already loaded
if isfield(handles,'spikes'), handles=wc_plot_spikes_and_ISI(handles); end
guidata(handles.mainfig, handles);

%% CLuster buttons

function rejectbuttons_Callback(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
% Reject cluster i
tag=get(source,'Tag');
i=str2double(tag(7:end));

if size(handles.WC.clus_per_temp,2)>=i %% temporary debug Flaff 20160617 ch3-2-su
handles.WC.clus_per_temp(:,i)=[];
end

%% Move the fix and fuse checkbox information along
if find(handles.fixed),
    handles.fixed(i:end-1)=handles.fixed(i+1:end);
    handles.fixed(end)=0;
    handles.fixed_classind={handles.fixed_classind{setdiff(1:length(handles.fixed_classind),i)}};
    handles.fixed_classind{end+1}=[];
    for j=i:length(handles.hfix)-1,
        set(handles.hfix(j),'Value',get(handles.hfix(j+1),'Value'));
        set(handles.hfuse(j),'Value',get(handles.hfuse(j+1),'Value'));
    end
    set(handles.hfix(end),'Value',0);
    set(handles.hfuse(end),'Value',0);
end

%% Either create unsorted cluster or add to existing unsorted cluster
if numel(handles.classind)>handles.ncl
    handles.classind{end}=[handles.classind{end} handles.classind{i}];
else
    handles.classind{end +1}=handles.classind{i};
end
handles.classind=handles.classind(setdiff(1:length(handles.classind),i));
handles.classind_unforced{i}=[];

handles.ncl=handles.ncl-1;
handles.rejected=i;
handles=wc_plot_spikes_and_ISI(handles);
handles=wc_plot_temperature(handles);
handles.rejected=1;

guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function fusebutton_Callback(source,~)
% Fuse selected clusters cluster
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
to_fuse=[];
for i=1:numel(handles.hfuse)
    if isobject(handles.hfuse(i))
        if get(handles.hfuse(i),'Value')
            to_fuse=[to_fuse, i];
        end
    end
end
to_fuse(to_fuse>numel(handles.classind))=[];
if numel(to_fuse)>1
    f=min(to_fuse);
    if size(handles.WC.clus_per_temp,2)>=to_fuse(end)
        clustemp=handles.WC.clus_per_temp(:,to_fuse);
        to_fuse=to_fuse(to_fuse~=f);
        handles.WC.clus_per_temp(:,to_fuse)=[];
    else
        clustemp=[handles.WC.clus_per_temp(:,to_fuse(1:end-1)) handles.WC.clus_per_temp(:,f)];
        to_fuse=to_fuse(to_fuse~=f);
        handles.WC.clus_per_temp(:,to_fuse(1:end-1))=[];
    end
    already_fused=0; %keep track of number of shifts
    M=numel(handles.classind{f});
    largest_cluster=f; 
    for i=to_fuse
        i=i-already_fused;
        n=numel(handles.classind{i});
        handles.classind{f}         =[handles.classind{f} handles.classind{i}];
        handles.classind            =handles.classind(1:length(handles.classind)~=i);
        handles.classind_unforced{f}=[handles.classind_unforced{f} handles.classind_unforced{i}];        
        handles.classind_unforced   =handles.classind_unforced(1:length(handles.classind_unforced)~=i);
        handles.fixed_classind=handles.fixed_classind(1:length(handles.fixed_classind)~=i);
        handles.fixed_classind{end+1}=[];
        %shift fix checkbox values
        for j=i:numel(handles.fixed)-1
            handles.fixed(j)=handles.fixed(j+1);
            set(handles.hfix(j),'Value',get(handles.hfix(j+1),'Value'));
        end
        % to keep cluster markers for the largest of fused clusters
        if n>M
            largest_cluster=i+already_fused;
            M=n;
        end
        already_fused=already_fused+1;
    end    
    handles.WC.clus_per_temp(:,f)=clustemp(:,[f to_fuse]==largest_cluster);
    n_fused=numel(to_fuse);
    handles.fixed(end-n_fused+1:end)=0;
    set(handles.hfix(end-n_fused+1:end),'Value',0);
    handles.ncl=handles.ncl-n_fused;
end

%reset fuse
for i=1:numel(handles.hfuse)
    set(handles.hfuse(i),'Value',0);
end
handles=wc_plot_spikes_and_ISI(handles);
handles=wc_plot_temperature(handles);
guidata(handles.mainfig, handles); % Update handles structure
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function fusebuttons_Callback(source,~)

function fixbuttons_Callback(source,~)
handles=guidata(get(source,'UserData'));
tag=get(source,'Tag');
i=str2double(tag(4:end));
switch get(source,'value');
    case 1, handles.fixed(i)=1;handles.fixed_classind{i}=handles.classind{i};
    case 0, handles.fixed(i)=0;handles.fixed_classind{i}=[];
end
guidata(handles.mainfig, handles); % Update handles structure

function restclusterbutton_Callback(source,~)
handles=guidata(get(source,'UserData'));
handles.ncl=numel(handles.classind);
handles.classind_unforced{end+1}=handles.classind_unforced{end}; %% end has to be the new cluster, end+1 has to be the rest from not sorted
handles.classind_unforced{end-1}=[];
handles.classind{end+1}=[]; 
handles=wc_plot_spikes_and_ISI(handles);
handles.WC.clus_per_temp=[handles.WC.clus_per_temp [2;2]];
guidata(handles.mainfig, handles);

%% broadband functions

function previous_window(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
window_length=handles.ts_time(2)-handles.ts_time(1);
handles.ts_time=handles.ts_time-window_length;
set_axes_limits(handles);
set(handles.start_t_textbox,'String',num2str(handles.ts_time(1)));
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function next_window(source,~)
set(handles.textStatus,'string',sprintf('Processing...'));
handles=guidata(get(source,'UserData'));
window_length=handles.ts_time(2)-handles.ts_time(1);
handles.ts_time=handles.ts_time+window_length;
set_axes_limits(handles);
set(handles.start_t_textbox,'String',num2str(handles.ts_time(1)));
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function apply_start_t(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
start_t=str2double(get(source,'String'));
window_length=handles.ts_time(2)-handles.ts_time(1);
handles.ts_time(1)=start_t;
handles.ts_time(2)=start_t+window_length;
set_axes_limits(handles);
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function apply_window_length(source,~)
    handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
    window_length=str2double(get(source,'String'));
    handles.ts_time(2)=handles.ts_time(1)+window_length;
    set_axes_limits(handles);
    guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

function set_axes_limits(handles)
set(handles.axesTS,'xlim',[min(handles.ts_time) max(handles.ts_time)]);
set(handles.axesTS,'ylim',[-handles.WC.thr(1) handles.WC.thr(1)]*3);

function plot_timecourse(source,~)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
ifplot=get(handles.plot_timecourse,'Value');
cla(handles.axesTS);
if ifplot
    filename= handles.filename;
    pathname=handles.pathname;
    load([pathname filesep 'concatenation_info.mat'],'blocksamplesperchannel','wheretofindwhat','whattofindwhere','channels_to_process','sr');
    us_idx=strfind(filename,'_');
    %n_file=str2double(filename(us_idx(2)+1:us_idx(3)-1));
    n_file=str2double(filename(us_idx(2)+1:end));
    blocks=whattofindwhere{handles.channel}{n_file};
    us_idx=strfind(pathname,filesep);
    BB_data=[];
    for b=1:numel(blocks)
        block=blocks(b);
        tens_fname=[pathname(1:us_idx(end-1)) 'WC_Block-' num2str(block) filesep 'datafilt_ch' sprintf('%03d.mat',handles.channel)];
        load(tens_fname);
        BB_data=[BB_data data];
        BB_block_samples(b)=numel(BB_data);
    end
    
    wc_bins=1/handles.WC.sr:1/handles.WC.sr:numel(BB_data)/handles.WC.sr;    
    set(handles.textStatus,'string',sprintf('Plotting %s',handles.filename));    
    hold(handles.axesTS,'on');
    plot(handles.axesTS,wc_bins, BB_data*handles.WC.transform_factor,'color','k','parent',handles.axesTS);
    plot(handles.axesTS,[wc_bins(1);wc_bins(end)], repmat(-handles.WC.thr(1),[2 1]),'color','r','parent',handles.axesTS);
    plot(handles.axesTS,[wc_bins(1);wc_bins(end)], repmat(handles.WC.thr(1),[2 1]),'color','r','parent',handles.axesTS);
    plot(handles.axesTS,[wc_bins(1);wc_bins(end)], repmat(-handles.WC.thr(2),[2 1]),'color','g','parent',handles.axesTS);
    plot(handles.axesTS,[wc_bins(1);wc_bins(end)], repmat(handles.WC.thr(2),[2 1]),'color','g','parent',handles.axesTS);
    handles.spike_indicators=plot(handles.axesTS,[handles.index/1000 handles.index/1000]',repmat([3;2]*-handles.WC.thr(1),1,numel(handles.index)),'k');    
    set(handles.axesTS,'xlim',[min(handles.ts_time) max(handles.ts_time)]);
    set(handles.axesTS,'ylim',[-handles.WC.thr(1) handles.WC.thr(1)]*3);
end
handles=wc_plot_spikes_and_ISI(handles); %% needs to be her so spike indicator colors are always updated as well
guidata(handles.mainfig, handles);
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));

%% main figure (and additional cluster) functions

function handles=create_mainfig(handles)
handles.isaGUI=1;
handles.mainfig=figure('Visible','on','Units','Normalized','Position',[0 0 1,0.9],...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Name','Wave_clus','NumberTitle','off',...
    'Tag','Mainfig');
handles.ts_time=[0 10]; % this needs to be set somehow
figure(handles.mainfig);

ncol=5;
nrow=4;
stepx=0.02; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;


%% bb part
handles.ts_time=[0 10];
stepy2=0.03;
y_pos_ts=1-(stepy+hight)*1.05;
uicontrol('Style','pushbutton','String','Previous Window', 'units','normalized','Position',[0.7,y_pos_ts,0.08,stepy2],'UserData',handles.mainfig,'Callback',@previous_window);
uicontrol('Style','pushbutton','String','Next Window',     'units','normalized','Position',[0.8,y_pos_ts,0.08,stepy2],'UserData',handles.mainfig,'Callback',@next_window);
uicontrol('Style','text','String','Start (s)','units','normalized','Position',[0.5,y_pos_ts,0.04,stepy2]);
handles.start_t_textbox       = uicontrol('Style','edit','units','normalized','String',num2str(handles.ts_time(1)),'Position',[0.55,y_pos_ts,0.04,stepy2],'UserData',handles.mainfig,'Callback',@apply_start_t);
uicontrol('Style','text','String','Window (s)','units','normalized','Position',[0.6,y_pos_ts,0.04,stepy2]);
handles.window_length_textbox = uicontrol('Style','edit','units','normalized','String',num2str(diff(handles.ts_time)),'Position',[0.65,y_pos_ts,0.04,stepy2],'UserData',handles.mainfig,'Callback',@apply_window_length);
handles.plot_timecourse = uicontrol('units','normalized','Style','checkbox','String','Plot timecourse','FontSize',12,'Tag','Plot_t',...
            'UserData',handles.mainfig,'Position',[0.4,y_pos_ts,0.08,stepy2],'Callback',{@plot_timecourse});
%% plot all spikes (for feature plots only so far)
handles.plot_all = uicontrol('units','normalized','Style','checkbox','String','All','FontSize',12,'Tag','Plot_t',...
            'UserData',handles.mainfig,'Position',[stepx+0.11 1-(stepy+hight)*( 1 )-stepy*0.9 0.02 0.03],'Callback',@plot_all_spikes);
        
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
    'Position',[stepx+(width+stepx)*1.5 1-stepy*0.9 0.35 0.035]);
handles.hsave=uicontrol('units','normalized','Style','pushbutton','String','Save','FontSize',12,...
    'Tag','save',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*3.5 1-stepy*0.9 0.05 0.035],...
    'Callback',{@savebutton_Callback});
handles.hnext=uicontrol('units','normalized','Style','pushbutton','String','Next','FontSize',12,...
    'Tag','Next',...
    'UserData',handles.mainfig,...
    'Position',[stepx+(width+stepx)*4 1-stepy*0.9 0.05 0.035],...
    'Callback',{@nextbutton_Callback});
handles.axesTS=axes('position',[stepx 1-(stepy+hight)*( 1 )+stepy*0.5 1-2*stepx hight*0.9],...
    'Tag','TS',...
    'UserData',handles.mainfig,...
    'ButtonDownFcn',{@copy_to_new_window});
handles.t=1;
handles.axesAllClusters=axes('position',[stepx 1-(stepy+hight)*( 2 ) width hight],...load
    'Tag','AllClusters','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.axesClust0=axes('position',[stepx+(width+stepx)*0 1-(stepy+hight)*( nrow ) width hight],...
    'Tag','Clust0','ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
handles.hclassify=uicontrol('units','normalized','Style','togglebutton','String','Class T','FontSize',12,...
    'Tag','Class T',...
    'UserData',handles.mainfig,...
    'Position',[stepx-0.01  1-(stepy+hight)*( nrow )-stepy*0.9 0.04 0.035],...
    'Callback',{@classifybutton_Callback});
handles.htempmatch=uicontrol('units','normalized','Style','togglebutton','String','Classify','FontSize',12,...
    'Tag','Classify',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.03 1-(stepy+hight)*( nrow )-stepy*0.9 0.03 0.035],...
    'Callback',{@classifybutton_Callback});
handles.htempmatch2=uicontrol('units','normalized','Style','togglebutton','String','Near','FontSize',12,...
    'Tag','Near',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.06 1-(stepy+hight)*( nrow )-stepy*0.9 0.03 0.035],...
    'Callback',{@classifybutton_Callback});
handles.htempmatch3=uicontrol('units','normalized','Style','togglebutton','String','Near T','FontSize',12,...
    'Tag','NearT',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.09 1-(stepy+hight)*( nrow )-stepy*0.9 0.03 0.035],...
    'Callback',{@classifybutton_Callback});
handles.hfuseclusters=uicontrol('units','normalized','Style','togglebutton','String','Fuse','FontSize',12,...
    'Tag','fuse',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.12 1-(stepy+hight)*( nrow )-stepy*0.9 0.03 0.035],...
    'Callback',{@fusebutton_Callback});
handles.hrest=uicontrol('units','normalized','Style','pushbutton','String','-->','FontSize',12,...
    'Tag','-->',...
    'UserData',handles.mainfig,...
    'Position',[stepx+0.15 1-(stepy+hight)*( nrow )-stepy*0.9 0.03 0.035],...
    'Callback',{@restclusterbutton_Callback});

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
    'Position',[stepx+0.14 1-(stepy+hight)*( 1 )-stepy*0.9 0.04 0.03],...
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
        x=stepx+(width+stepx)*c;
        y= 1-(stepy+hight)*r-stepy*0.45;
        handles.spikeaxes(i)=axes('position',[x y+stepy*0.45 width hight],...
            'Tag',sprintf('Clust%d',i),'ButtonDownFcn',{@copy_to_new_window},'NextPlot','add');
        handles.hfuse(i) = uicontrol('units','normalized','Style','checkbox','String','Fuse','FontSize',12,...
            'Tag',sprintf('Fuse%d',i),'UserData',handles.mainfig,'Position',[x+0.03 y 0.035 0.02],'Callback',{@fusebuttons_Callback});
        handles.hfix(i) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
            'Tag',sprintf('Fix%d',i),'UserData',handles.mainfig,'Position',[x y 0.025 0.02],'Callback',{@fixbuttons_Callback});
        handles.hreject(i) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
            'Tag',sprintf('Reject%d',i),'UserData',handles.mainfig,'Position',[x+0.08 y 0.01 0.02],'Callback',{@rejectbuttons_Callback});
        handles.hdetails(i) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
            'Tag',sprintf('Details%d',i),'UserData',handles.mainfig,'Position',[x+0.11 y 0.06 0.02],'Callback',{@detailsbuttons_Callback});
        handles.hclustergroup{i}=[handles.spikeaxes(i) handles.hfuse(i) handles.hfix(i) handles.hreject(i) handles.hdetails(i)];
        set(handles.hclustergroup{i},'Visible','Off');
    end
end

function handles=create_supplfig(handles)
handles.hsuppl = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig','CloseRequestFcn','close_or_make_invisible');
nrow=4;ncol=5;
stepx=0.02; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;
j=12;
for r=1:nrow
    for c=1:ncol,
        j=j+1;
        x=stepx+(width+stepx)*(c-1);
        y= 1-(stepy+hight)*r-stepy*0.45;
        
        handles.spikeaxes(j)=axes('position',[x y+stepy*0.75 width hight],...
            'Tag',sprintf('Clust%d',j));
        handles.hfuse(j) = uicontrol('units','normalized','Style','checkbox','String','Fuse','FontSize',12,...
            'Tag',sprintf('Fuse%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.03 y 0.035 0.02],...
            'Callback',{@fusebuttons_Callback});
        handles.hfix(j) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
            'Tag',sprintf('Fix%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x y  0.025 0.02],...
            'Callback',{@fixbuttons_Callback});
        handles.hreject(j) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
            'Tag',sprintf('Reject%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.08 y 0.01 0.02],...
            'Callback',{@rejectbuttons_Callback});
        handles.hdetails(j) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
            'Tag',sprintf('Details%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.11 y 0.06 0.02],...
            'Callback',{@detailsbuttons_Callback});
        handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfuse(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j)];
        set(handles.hclustergroup{j},'Visible','Off');
    end
end

function handles=create_supplfig2(handles)
handles.hsuppl2 = figure('Visible','Off','Units','Normalized','Position',[0 0 1 0.9],...
    'Name','More clusters','NumberTitle','off',...
    'Paperunits','points','Paperorientation','portrait','PaperPosition',[0 0 1920 1080],...
    'Tag','Supplfig','CloseRequestFcn','close_or_make_invisible');
nrow=4;ncol=5;
stepx=0.02; width=(1-(ncol+1)*stepx)/ncol;
stepy=0.05; hight=(1-(nrow+1)*stepy)/nrow;
j=32;
for r=1:nrow
    for c=1:ncol,
        j=j+1;
        x=stepx+(width+stepx)*(c-1);
        y= 1-(stepy+hight)*r-stepy*0.45;
        
        handles.spikeaxes(j)=axes('position',[x y+stepy*0.75 width hight],...
            'Tag',sprintf('Clust%d',j));
        handles.hfuse(j) = uicontrol('units','normalized','Style','checkbox','String','Fuse','FontSize',12,...
            'Tag',sprintf('Fuse%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.03 y 0.035 0.02],...
            'Callback',{@fusebuttons_Callback});
        handles.hfix(j) = uicontrol('units','normalized','Style','checkbox','String','Fix','FontSize',12,...
            'Tag',sprintf('Fix%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x y  0.025 0.02],...
            'Callback',{@fixbuttons_Callback});
        handles.hreject(j) = uicontrol('units','normalized','Style','pushbutton','String','X','FontSize',14,'ForeGroundColor','r',...
            'Tag',sprintf('Reject%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.08 y 0.01 0.02],...
            'Callback',{@rejectbuttons_Callback});
        handles.hdetails(j) = uicontrol('units','normalized','Style','pushbutton','String','Details','FontSize',12,...
            'Tag',sprintf('Details%d',j),...
            'UserData',handles.mainfig,...
            'Position',[x+0.11 y 0.06 0.02],...
            'Callback',{@detailsbuttons_Callback});
        handles.hclustergroup{j}=[handles.spikeaxes(j) handles.hfuse(j) handles.hfix(j) handles.hreject(j) handles.hdetails(j)];
        set(handles.hclustergroup{j},'Visible','Off');
    end
end

function close_or_make_invisible(source, ~) %not used??
handles=guidata(get(source,'UserData'));
if isvalid(handles.mainfig)
    set(source,'visible','off');
else
    close(source);
end
guidata(handles.mainfig, handles); 

function plot_all_spikes(source, event)
handles=guidata(get(source,'UserData'));
set(handles.textStatus,'string',sprintf('Processing...'));
all=get(source,'Value');
if all
handles.const_MAX_SPIKES_TO_PLOT=handles.nspk;
else
handles.const_MAX_SPIKES_TO_PLOT=1000;
end
handles=wc_plot_spikes_and_ISI(handles);
guidata(handles.mainfig, handles); 
set(handles.textStatus,'string',sprintf('%s',[handles.pathname(numel(handles.datafolder)+1:end) handles.filename]));


function copy_to_new_window(source,~) %kinda useless too? --> maybe used for details??
figure(100);
cla;
copyobj(get(source,'children'),gca);
