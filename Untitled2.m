
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
    
    n_classes=0;
    for i=1:handles.WC.max_nrclasses,
        classind{i}=[];
    end
    temp_start=1;
    handles.WC.temp=[];
    handles.WC.n_clus_per_temp={};
    while n_classes < handles.WC.max_nrclasses && temp_start<size(tree,1)% max clusters ( leave 1 for unclustered) not reached yet
        [temp] = wc_find_temperature(tree(temp_start:end,:),handles); %% need to reduce tree?
        temp=temp+temp_start-1;
        
        %DEFINE CLUSTERS for specific temperature using min_clus variable
        %             classind=cell(1,handles.WC.max_nrclasses);
        %             t_remove=true(1,handles.WC.max_nrclasses);
        n_classes=sum(~cellfun(@isempty,classind));
        clusters_for_this_temp=[];
        for i=2:handles.WC.max_nrclasses,
            t=classtemp{temp,i};
            t = setdiff(t,[classind{:}]);
            if length(t)>min_clus && n_classes < handles.WC.max_nrclasses
                n_classes=n_classes+1;
                classind{n_classes}=t;
                clusters_for_this_temp=[clusters_for_this_temp i];
            end
        end
        if temp==temp_start % didnt find appropriate temperature
            break
        end
        if clusters_for_this_temp
            handles.WC.temp=[handles.WC.temp temp]; %% store all used temperatures
            handles.WC.n_clus_per_temp=[handles.WC.n_clus_per_temp {clusters_for_this_temp}]; %% store all used temperatures
        end
        temp_start=temp;
    end
    %         if isempty(handles.WC.n_clus_per_temp)
    %             handles.WC.n_clus_per_temp=0;
    %         end
    %% add biggest cluster of last iteration
    t=classtemp{temp,1};
    t = setdiff(t,[classind{:}]);
    if n_classes+1<handles.WC.max_nrclasses;
        classind{n_classes+1}=t;
        NN=numel(handles.WC.n_clus_per_temp);
        if NN>0
            handles.WC.n_clus_per_temp{NN}= [handles.WC.n_clus_per_temp{NN}, 1];
        else
            handles.WC.n_clus_per_temp{1}= 1;
        end
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
    
    %% Forcing as default
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
            %eval(sprintf('print -djpeg fig2print_%s-%d',handles.bname,channelfile));
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
    %    memory
end