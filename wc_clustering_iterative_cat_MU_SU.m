function wc_clustering_iterative_cat_MU_SU(handles)
handles.bname='Ch';
ifplot=1;
print2file =1;                              %for saving printouts. ??

if ispc, handles.system = 'windows'; elseif ismac, handles.system = 'MACI64'; else if isunix, handles.system = 'linux'; else error('TestClust:UnknownSystem','Unknown system');
    end
end
sr = handles.WC.sr;

clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0 0 0.5;0 0 0];
set(0,'DefaultAxesColorOrder',clus_colors)

switch handles.WC.threshold
    case 'pos'
        thresholds={[handles.current_threshold_step '_pos']};
    case 'neg'
        thresholds={[handles.current_threshold_step '_neg']};
    case 'both'
        thresholds={[handles.current_threshold_step '_neg'],[handles.current_threshold_step '_pos']};
end
redoCount = ones(numel(thresholds),1);


for k =  1 : numel(thresholds)
    
    %    memory
    tic
    channelfile=[sprintf('%03d',handles.current_channel) '_' num2str(handles.current_channel_file)];
    filename=['dataspikes_ch' channelfile '_' thresholds{k} '.mat'];
    handles.channel=handles.current_channel;
    handles.thrsign=thresholds{k};
    
    disp(filename)
    load([handles.WC_concatenation_folder filename])
    
    if isempty(spikes), fprintf('No spikes detected in file %s\n',channelfile); continue; end
    nspk = size(spikes,1);
    if nspk<15, fprintf('Less than 15 spikes detected in file %s\n',channelfile); continue; end
    
    fprintf('%d spikes detected\n',nspk);
    
    handles.fname = [handles.WC_concatenation_folder 'data_' channelfile '_' thresholds{k}];   %filename for interaction with SPC
    handles.fname = ['data_' channelfile '_' thresholds{k}];   %filename for interaction with SPC
    
    min_clus = max(handles.WC.min_clus_abs,handles.WC.min_clus_rel*min(nspk,handles.WC.max_spikes2cluster));
    handles.WC.min_clus = min_clus;
    
    fprintf('Feature detection...\n');
        yyy=any(strfind(filename,'dataspikes_ch001_2_SU_neg'));
    %CALCULATES INPUTS TO THE CLUSTERING ALGORITHM
    [features,feature_names,inputs] = wc_feature_selection(spikes,index,handles);
    handles.WC.inputs = inputs;
    save([handles.WC_concatenation_folder filename],'features','feature_names','-append');
    
    redo = 1;
    
    if ~isempty(features)
        
        while redo && redoCount(k)<100
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
                redoCount(k) = redoCount(k) + 1;
            end
        end
        classtemp=cell(size(clu,1),handles.WC.max_nrclasses);
        for te=1:size(clu,1) %temperature
            for j=1:handles.WC.max_nrclasses,
                classtemp{te,j}=inds_to_cluster((clu(te,3:end)==j-1));
            end;
        end
        fprintf('Adding SPC output information into results file...\n');
        save([handles.WC_concatenation_folder filename],'tree','classtemp','-append');
        
        
        %% automatic temperature selection
        for i=1:handles.WC.max_nrclasses,
            classind{i}=[];
        end
        
        min_clus=handles.WC.min_clus;
        
        reduced_tree=tree(2:end,5:end);
        diff_matrix=diff(tree(:,5:end));
        diff_matrix_ind=diff_matrix>min_clus & [diff_matrix(2:end,:); -1*ones(1,size(diff_matrix,2))]<0 & tree(2:end,5:end)>min_clus;
        [ts, clus]=ind2sub(size(diff_matrix_ind),find(diff_matrix_ind));
        clus_sizes=reduced_tree(sub2ind(size(reduced_tree),ts,clus));
        [~, sorted_idx]=sort(clus_sizes);
        sorted_idx=sorted_idx(numel(sorted_idx)-min(handles.WC.max_nrclasses-1,numel(sorted_idx))+1:end);
        
        handles.WC.clus_per_temp=[ts(sorted_idx)'+1; clus(sorted_idx)'];
        valid=[];
        n_classes=0;
        for i=1:size(handles.WC.clus_per_temp,2)
            t=classtemp{ts(i)+1,clus(i)};
            t = setdiff(t,[classind{:}]);
            if numel(t)>0
            n_classes=n_classes+1;
            classind{n_classes}=t;
            valid=[valid; i];
            end
        end
        handles.WC.clus_per_temp=handles.WC.clus_per_temp(:,valid);
        
        % add biggest cluster of last iteration
        if ~isempty(ts)
            temp_t=max(ts)+1;
        else
            temp_t=1;
        end
        t=classtemp{temp_t,1};
        t = setdiff(t,[classind{:}]);
        if n_classes+1<handles.WC.max_nrclasses;
            classind{n_classes+1}=t;            
            handles.WC.clus_per_temp=[handles.WC.clus_per_temp [temp_t; 1]];
        end
        classind(cellfun(@isempty,classind))=[];
        
        %zero cluster, includes all unclustered spikes
        classind{end+1}=setdiff(1:nspk, [classind{:}]);        
        
        %% classify rest as default
        handles.spikes=spikes;
        handles.classind=classind;
        handles.features=features;
        handles.feature_names=feature_names;
        handles.ncl=max(length(classind)-1);
        handles.nspk=nspk;
        clear features feature_names
        handles=wc_classifyrest(handles);
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
            handles.mainfig=[];
            handles=wc_plot_features_vs_time(handles);
            print(handles.htimecourse,[handles.bname '-' channelfile '_' thresholds{k} '-timecourse'],'-djpeg');
            close (handles.htimecourse)            
            
            %% export features figure
            handles.const_MAX_SPIKES_TO_PLOT=1000; %to prevent large plottings
            handles.index=index;
            handles.nfeatures=size(handles.features,2);
            set(0,'DefaultAxesColorOrder',clus_colors);
            handles.colors= clus_colors;
            handles.mainfig=[];
            handles=wc_plot_features_vs_features(handles);
            print(handles.hfeatures,[handles.bname '-' channelfile '_' thresholds{k} '-features'],'-djpeg');
            close (handles.hfeatures)
            
            %% plotting WC main figure
            handles.rejected=1; %%??
            handles.plotted=[];
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
                print(handles.mainfig,[handles.bname '-' channelfile '_' thresholds{k} '-clusters'],'-djpeg');
                close(handles.mainfig);
            end
        end
        %spikes file
        clear time
        clear index_ spikes_
        clear cluster_class
        clear clu tree classtemp
        clear classind cluster_class
    end
    disp([' clustering ' handles.WC_concatenation_folder filename ' took ' num2str(round(toc)) ' seconds']);        
end
clear handles
