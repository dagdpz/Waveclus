function Do_clustering4_redo(handles)

ifplot=0;
print2file =1;                              %for saving printouts.
%print2file =0;                              %for printing printouts.

if ispc, handles.system = 'windows'; elseif ismac, handles.system = 'MACI64'; else if isunix, handles.system = 'linux'; else error('TestClust:UnknownSystem','Unknown system');
    end
end

handles.par.features = 'wavpcarawderiv';    %choice of spike features: wav: wavelet decomposition; pca: principle component analyses; raw: raw waveforms; deriv: first derivative of the raw waveforms
handles.par.features = 'wavpcaraw';
handles.par.wavelet='haar';                 %choice of wavelet family for wavelet features
handles.par.exclusioncrit = 'thr';          % thr; number
handles.par.exclusionthr = 0.8;               %def R^2 = 0.80
handles.par.maxinputs = 15;   %17, 15              %number of inputs to the clustering def. 11
handles.par.scales = 4;                     %scales for wavelet decomposition

% handles.bname='datahighpassch';

%SPC parameters
handles.par.num_temp = 18;                  %number of temperatures; def 25
handles.par.mintemp = 0;                    %minimum temperature
handles.par.maxtemp = 0.18;                 %maximum temperature def 0.25
handles.par.tempstep = 0.01;                %temperature step
% handles.stab = 0.8;                       %stability condition for selecting the temperature
handles.par.SWCycles = 1000;  % def. 100     %number of montecarlo iterations
handles.par.KNearNeighb = 11;               %number of nearest neighbors

handles.par.temp_plot = 'log';              % temperature plot in log scale

handles.par.max_spikes2plot = 1000;         %maximum number of spikes to plot.

handles.par.min_clus_abs = 10;
handles.par.min_clus_rel = 0.005;          %Default: 0.005% alternative: 0.0035
handles.par.max_nrclasses = 13;
handles.par.max_nrclasses2plot = 8;

handles.par.chunk=5;                        %length of pieces into which file has to be splitted
handles.par.max_spikes2cluster = 30000;       %maximum number of spikes to cluster, if more take only this amount of randomly chosen spikes, others are set into cluster 0

sr = handles.par.sr;

clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5];
set(0,'DefaultAxesColorOrder',clus_colors)

channels = dir('dataspikes_ch*.mat');
channels = {channels.name};

redoCount = ones(length(channels),1);

for k =  1 : length(channels)
    
    disp('Resetting initial values')
    handles.par.SWCycles = 1000;
    handles.par.maxinputs = 15;
    handles.par.max_spikes2cluster = 30000;
    
    
    %    memory
    tic
    channel=channels{k}(1:end-4);
    handles.channel=str2double(channel(strfind(channel,'ch')+2:strfind(channel,'ch')+4));
    handles.thrsign=channel(strfind(channel,'ch')+6 : end-4);
    
    disp(channel)
    load(channel)
    
    if isempty(spikes), fprintf('No spikes detected in file %s\n',filename); continue; end
    nspk = size(spikes,1);
    if nspk<15, fprintf('Less than 15 spikes detected in file %s\n',filename); continue; end
    
    fprintf('%d spikes detected\n',nspk);
    
    handles.fname = ['data_' channel];   %filename for interaction with SPC
    
    min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
    handles.par.min_clus = min_clus;
    
    fprintf('Feature detection...\n');
    %CALCULATES INPUTS TO THE CLUSTERING ALGORITHM
    [features,feature_names,inputs] = wave_features5B(spikes,handles);
    handles.par.inputs = inputs;
    save(channel,'features','feature_names','-append');
    
    redo = 1;
    
    if ~isempty(features)
        
        while redo
            try
                %INTERACTION WITH SPC
                fprintf('SPC...\n');
                if nspk>handles.par.max_spikes2cluster, %take random handles.par.max_spikes2cluster spikes
                    t=randperm(nspk);
                    inds_to_cluster=t(1:handles.par.max_spikes2cluster);
                    features1=features(inds_to_cluster,:);
                else
                    features1=features;
                    inds_to_cluster=1:nspk;
                end
                
                save(handles.fname,'features1','-ascii');
                
                [clu, tree] = run_cluster(handles);
                
                redo = 0;
            catch
                disp('Again')
                redo = 1;
                redoCount(k) = redoCount(k) + 1;
                
                if handles.par.SWCycles > 200
                    handles.par.SWCycles = handles.par.SWCycles - 200;
                    disp(['New SWCycles: ' num2str(handles.par.SWCycles)])
                end
                if handles.par.maxinputs > 11
                    handles.par.maxinputs = handles.par.maxinputs - 2;
                    disp(['New maxinputs: ' num2str(handles.par.maxinputs)])
                end
                if handles.par.max_spikes2cluster > 20000
                    handles.par.max_spikes2cluster = handles.par.max_spikes2cluster - 10000;
                    disp(['New max_spikes2cluster: ' num2str(handles.par.max_spikes2cluster)])
                end
                disp('Recalculating Features')
                [features,feature_names,inputs] = wave_features5B(spikes,handles);
                handles.par.inputs = inputs;
                save(channel,'features','feature_names','-append');
            end
        end
        clear features feature_names
        
        classtemp=cell(size(clu,1),handles.par.max_nrclasses);
        for te=1:size(clu,1),
            %       for j=1:handles.par.max_nrclasses, classtemp{te,j}=inds_to_cluster(find(clu(te,3:end)==j-1)); end;
            for j=1:handles.par.max_nrclasses, classtemp{te,j}=inds_to_cluster((clu(te,3:end)==j-1)); end;
        end
        fprintf('Adding SPC output information into results file...\n');
        %    save(output_filename,'tree','classtemp','-append');
        %    save(output_filename,'spikes','-append');
        save(channel,'tree','classtemp','-append');
        
        [temp] = find_temp_new(tree,handles);
        handles.par.temp=temp;
        
        %DEFINE CLUSTERS for spesific temperature using min_clus variable
        classind=cell(1,handles.par.max_nrclasses);
        for i=1:handles.par.max_nrclasses,
            t=classtemp{temp,i};
            classind{i}=[];
            if length(t)>min_clus, classind{i}=t; end
        end
        %zero cluster, alo includes all unclustered spikes
        classind{end+1}=setdiff(1:nspk, [classind{:}]);
        
        %prepare to save
        cluster_class=zeros(nspk,2);
        cluster_class(:,2)= index(:);
        for i=1:length(classind)-1, %minus cluster zero
            cluster_class(classind{i},1)=i;
        end
        par=handles.par;
        %    save(output_filename,'cluster_class','par','-append');
        save(channel,'cluster_class','par','-append','-v6');
        
        if ifplot,
            %plotting
            time=[-handles.par.w_pre+1:1:handles.par.w_post]/sr*1000;%timie in ms
            figure
            set(gcf,'PaperOrientation','Landscape','PaperUnits','normalized','PaperPosition',[0.01 0.01 0.98 0.98])
            clf
            ncol=5;
            temp_plot=6;
            if length(classind)>4, nrow=4;sp_plot=[2 3 4 5 11 12 13 14];
            else nrow=2;sp_plot=[2 3 4]; end
            
            isi_plot=sp_plot+ncol;
            
            %temperature plot
            subplot(nrow,ncol,temp_plot);
            semilogy(1:handles.par.num_temp,tree(1:handles.par.num_temp,5:end));
            line([temp temp],[1 max(max(tree(:,5:end)))*1.1],'linestyle',':','color','k');
            line([1 handles.par.num_temp],[min_clus min_clus],'linestyle',':','color','k');
            ylim([1 max(max(tree(:,5:end)))*1.1])
            xlim([0.5 handles.par.num_temp+0.5]);
            
            %all spikes superimposed
            if length(classind)>1,
                mn=min(min(spikes([classind{1:end-1}],:)));
                mx=max(max(spikes([classind{1:end-1}],:)));
            else %only cluster 0
                mn=min(min(spikes([classind{1}],:)));
                mx=max(max(spikes([classind{1}],:)));
            end
            
            subplot(nrow,ncol,1);cla; hold on;
            if sum(cluster_class(:,1)==0), plot(time, spikes(cluster_class(:,1)==0,:),'color','k'); end
            for i=1:min(length(classind)-1, handles.par.max_nrclasses2plot),
                if ~isempty(classind{i}), plot(time, spikes(classind{i},:),'color',clus_colors(i,1:3)); end
            end
            xlim([min(time) max(time)]);
            ylim([mn mx]);
            title(sprintf('%s-%d',handles.bname,channel),'Fontsize',12,'interpreter','none')
            
            %individual clusters
            for i=1:min(length(classind)-1, handles.par.max_nrclasses2plot),%minus cluster 0
                subplot(nrow,ncol,sp_plot(i));hold on
                sp=spikes(classind{i},:);
                plot(time, sp,'color',clus_colors(i,1:3));
                m=mean(sp);
                s=std(sp);
                plot(time, m,'color','k','linewidth',2);
                plot(time, m+s,'color','k','linewidth',0.5);
                plot(time, m-s,'color','k','linewidth',0.5);
                xlim([min(time) max(time)]);
                ylim([mn mx]);
                title(sprintf('Cluster %d, #%d',i,length(classind{i})));
                
                %isi
                subplot(nrow,ncol,isi_plot(i));
                isi=diff(index(classind{i}));
                edges=0:1:100;
                if isempty(isi), isi=0; end
                [N]=histc(isi,edges);
                h=bar(edges,N,'histc');
                set(h,'facecolor',clus_colors(i,1:3),'edgecolor',clus_colors(i,1:3),'linewidth',0.01);
                xlim([0 100]);
                title(sprintf('%d in <2ms',sum(N(1:2))));
                xlabel('ISI(ms)');
            end
            
            %cluster zero
            if sum(cluster_class(:,1)==0),
                subplot(nrow,ncol,sp_plot(end)+1);hold on
                sp=spikes(cluster_class(:,1)==0,:);
                plot(time, sp,'color','k');
                m=mean(sp);
                s=std(sp);
                plot(time, m,'color','c','linewidth',2);
                plot(time, m+s,'color','c','linewidth',0.5);
                plot(time, m-s,'color','c','linewidth',0.5);
                title(sprintf('Cluster 0, #%d',sum(cluster_class(:,1)==0)));
                xlim([min(time) max(time)]);
                %isi
                subplot(nrow,ncol,isi_plot(end)+1);
                isi=diff(index(cluster_class(:,1)==0));
                edges=0:1:100;
                if isempty(isi), isi=0; end
                [N]=histc(isi,edges);
                h=bar(edges,N,'histc');
                set(h,'facecolor','k','edgecolor','k','linewidth',0.01);
                xlim([0 100]);
                title(sprintf('%d in <2ms',sum(N(1:2))));
                xlabel('ISI(ms)');
            end
            toc
            if print2file==0;
                print
            else
                eval(sprintf('print -djpeg fig2print_%s-%d',handles.bname,channel));
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
    
end
clear handles
try
    figure;hist(redoCount,0:1:max(redoCount))
catch
    disp('nothing to plot')
end
