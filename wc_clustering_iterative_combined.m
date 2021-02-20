function wc_clustering_iterative_combined(features,feature_names,inputs,handles)
handles.bname='Ch';
ifplot=1;
print2file =1;                              %for saving printouts. ??

if ispc, handles.system = 'windows'; elseif ismac, handles.system = 'MACI64'; else if isunix, handles.system = 'linux'; else error('TestClust:UnknownSystem','Unknown system');
    end
end
sr = handles.WC.sr;

clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0 0 0.5;0 0 0];
set(0,'DefaultAxesColorOrder',clus_colors)

redoCount = 1;

tic
channelfile=[sprintf('%03d',handles.current_channel) '_' num2str(handles.current_channel_file)];
filename=['dataspikes_ch' channelfile '.mat'];
handles.channel=handles.current_channel;

disp(filename)
load([handles.WC_concatenation_folder filename])

if isempty(spikes), fprintf('No spikes detected in file %s\n',channelfile); return; end
nspk = size(spikes,1);
if nspk<15, fprintf('Less than 15 spikes detected in file %s\n',channelfile); return; end

fprintf('%d spikes detected\n',nspk);

handles.fname = [handles.WC_concatenation_folder 'data_' channelfile];   %filename for interaction with SPC
handles.fname = ['data_' channelfile];   %filename for interaction with SPC

min_clus = max(handles.WC.min_clus_abs,handles.WC.min_clus_rel*min(nspk,handles.WC.max_spikes2cluster));
handles.WC.min_clus = min_clus;

handles.WC.inputs = inputs;
save([handles.WC_concatenation_folder filename],'features','feature_names','-append');

redo = 1;

if ~isempty(features)
    while redo && redoCount<100
        try
            %INTERACTION WITH SPC
            fprintf('SPC...\n');
            if nspk>handles.WC.max_spikes2cluster, %take random handles.WC.max_spikes2cluster spikes
                t=randperm(nspk);
                inds_to_cluster=t(1:handles.WC.max_spikes2cluster);
                features1=features(inds_to_cluster,:);
            else
                features1=features;
                inds_to_cluster=1:nspk;
            end
            save(handles.fname,'features1','-ascii');
            [clu, tree] = wc_run_cluster(handles);
            redo = 0;
        catch err
            disp('Again')
            redo = 1;
            redoCount = redoCount + 1;
        end
    end
    classtemp=cell(size(clu,1),handles.WC.max_nrclasses);
    for te=1:size(clu,1) %temperature
        for j=1:handles.WC.max_nrclasses+1,
            classtemp{te,j}=inds_to_cluster((clu(te,3:end)==j-1));
        end;
    end
    fprintf('Adding SPC output information into results file...\n');
    save([handles.WC_concatenation_folder filename],'tree','classtemp','-append');
    
    for i=1:handles.WC.max_nrclasses,
        classind{i}=[];
    end
    [handles,classind]=wc_automatic_temperature_selection(handles,classind,classtemp,tree);
    
    
    
    %% classify rest as default
    handles.spikes=spikes;
    handles.classind=classind;
    handles.features=features;
    handles.feature_names=feature_names;
    handles.ncl=max(length(classind)-1);
    handles.nspk=nspk;
    handles.fixed=zeros(1,handles.ncl);
    clear features feature_names
    handles=wc_classifyrest(handles,1);
    classind=handles.classind;
    
    %prepare to save
    cluster_class=zeros(nspk,2);
    cluster_class(:,2)= index(:);
    for i=1:length(classind)-1, %minus cluster zero
        cluster_class(classind{i},1)=i;
    end
    par=handles.WC; %% saving as 'par'
    save([handles.WC_concatenation_folder filename],'cluster_class','par','-append','-v6');
    
    %% PLOTTING
    if ifplot,
        handles.isaGUI=0;
        %% export timecourse figure
        handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
        handles.index=index;
        handles.nfeatures=size(handles.features,2);
        set(0,'DefaultAxesColorOrder',clus_colors);
        handles.colors= clus_colors;
        handles.mainfig=99;
        handles=wc_plot_features_vs_time(handles);
        print(handles.htimecourse,[handles.bname '-' channelfile '-timecourse'],'-djpeg');
        close (handles.htimecourse)
        
        %% export features figure
        handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
        handles.index=index;
        handles.nfeatures=size(handles.features,2);
        set(0,'DefaultAxesColorOrder',clus_colors);
        handles.colors= clus_colors;
        handles.mainfig=99;
        handles=wc_plot_features_vs_features(handles);
        print(handles.hfeatures,[handles.bname '-' channelfile '-features'],'-djpeg');
        close (handles.hfeatures)
        
        %% plotting WC main figure
        handles.rejected=1; %%??
        handles.plotted=[];
        handles.isaGUI=0;
        handles.sp_time=1:(handles.WC.w_pre+handles.WC.w_post);
        handles=wc_create_mainfig(handles);
        handles=wc_plot_spikes_and_ISI(handles);
        
        handles.tree=tree;
        handles.min_clus=min_clus;
        handles=wc_plot_temperature(handles);
        
        set(gcf,'PaperUnits','normalized','PaperPosition',[0.01 0.01 0.98 0.98])
        if print2file==0;
            print
        else
            print(handles.mainfig,[handles.bname '-' channelfile '-clusters'],'-djpeg');
            close(handles.mainfig);
        end
    end
    %spikes file
%     clear time
%     clear index_ spikes_
%     clear cluster_class
%     clear clu tree classtemp
%     clear classind cluster_class
end
disp([' clustering ' handles.WC_concatenation_folder filename ' took ' num2str(round(toc)) ' seconds']);
clear handles
