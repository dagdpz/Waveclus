function Do_clustering_ION_pp(bname,channels,varargin)
% function Do_clustering_ION_long(bname,channels,filterset,posneg,prepost,whichwav,nfeatures)

var_names={'sr','filterset','posneg','prepost','features','whichwav','nfeatures','max_spikes2cluster','max_minutes2cluster','invert'};
defaults={25000,0,'neg',[10 22],'wavpca','haar',11,50000,Inf,0};
% gain=5*1e3/2048/1000;
p=parse_parameters(varargin,var_names,defaults);


ifplot=0;
print2file =1;                              %for saving printouts.
%print2file =0;                              %for printing printouts.

if ispc, handles.system = 'windows'; elseif ismac, handles.system = 'MACI64'; else if isunix, handles.system = 'linux'; else error('TestClust:UnknownSystem','Unknown system');
   end
end
% if exist('filterset','var'), handles.filterset=filterset; else handles.filterset=0; end
handles.filterset=p.filterset;

% if exist('posneg','var'), handles.threshold=posneg; else handles.threshold='neg'; end
handles.threshold=p.posneg;

% if exist('prepost','var'),
%    handles.par.w_pre=prepost(1);                       %number of pre-event data points stored
%    handles.par.w_post=prepost(2);                      %number of post-event data points stored
%    t=sum(prepost(1:2));
%    if round(log2(t))~=log2(t), error('TestClust:TimeBase','Spike time base has to be power of 2');end
% else
%    handles.par.w_pre=10;                       %number of pre-event data points stored
%    handles.par.w_post=22;                      %number of post-event data points stored
% end
handles.par.w_pre=p.prepost(1);                       %number of pre-event data points stored
handles.par.w_post=p.prepost(2);                      %number of post-event data points stored
t=sum(p.prepost(1:2));
if round(log2(t))~=log2(t), error('TestClust:TimeBase','Spike time base has to be power of 2');end

% handles.par.features = 'wavpca';            %choice of spike features
handles.par.features = p.features;
% if exist('whichwav','var'), handles.par.wavelet=whichwav;
% else handles.par.wavelet='haar'; end                 %choice of wavelet family for wavelet features
handles.par.wavelet=p.whichwav;

% if exist('nfeatures','var'), handles.par.inputs=nfeatures; 
% else
%     handles.par.inputs = 11;                    %number of inputs to the clustering
% end
handles.par.inputs=p.nfeatures;

handles.bname=sprintf('%s_w%d',bname,handles.par.w_pre+handles.par.w_post);

handles.par.stdmin = 5;                     %minimum threshold
handles.par.stdmax = 50;                    %maximum threshold
handles.par.interpolation = 'y';            %interpolation for alignment
handles.par.int_factor = 2;                 %interpolation factor

handles.par.scales = 4;                     %scales for wavelet decomposition

%SPC parameters
handles.par.num_temp = 25;                  %number of temperatures
handles.par.mintemp = 0;                    %minimum temperature
handles.par.maxtemp = 0.25;                 %maximum temperature
handles.par.tempstep = 0.01;                %temperature step
% handles.stab = 0.8;                       %stability condition for selecting the temperature
handles.par.SWCycles = 100;                 %number of montecarlo iterations
handles.par.KNearNeighb = 11;               %number of nearest neighbors

handles.par.temp_plot = 'log';              % temperature plot in log scale

handles.par.max_spikes2plot = 1000;         %maximum number of spikes to plot.

handles.par.min_clus_abs=10;
handles.par.min_clus_rel=0.005;%0.5%
handles.par.max_nrclasses=13;
handles.par.max_nrclasses2plot=8;

handles.par.transform_factor=(5*1e+6/(2048*20*1e+3));
handles.par.chunk=5;                        %length of pieces into which file has to be splitted
% handles.par.max_spikes2cluster=50000;       %maximum number of spikes to cluster, if more take only this amount of randomly chosen spikes, others are set into cluster 0
handles.par.max_spikes2cluster=p.max_spikes2cluster;
handles.par.sr = p.sr;
sr=p.sr;

clus_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.5 0 0; 0 0.5 0; 0 0 0.5;];
set(0,'DefaultAxesColorOrder',clus_colors)

switch handles.filterset
   case {-1,'1_5000'}, %no additional filters at all
      fmin=15000;
   case {0,'none'}, %no additional filters at all filters at all
      fmin = 0;
   case {1,'unfilt1000',1000},
      fmin = 1000;
   case {2,'unfilt300',300},
      fmin = 300;
   case {3,'unfilt300filt1000'},
      fmin=3001000;
   case {4,'unfilt300filt1000unfilt1000'},
      fmin=311000;
   case {5,'filt300'}, 
      fmin=3000;
   case {6,'filt1000'}, 
      fmin=10000;
end

for k=1:length(channels),
    %    memory
   tic
   channel=channels(k);
   handles.channel=channel;
   filename=sprintf('%s-%d.sp',bname,channel);
   handles.filename=filename;
   if ~exist(filename,'file'),  fprintf('File %s does not exist\n',filename); continue; end

   % Extract Spikes DATA
   %split file into chuncks of minutes
   f=fopen(filename,'r');
   fbeg=0;
   fseek(f,fbeg,'eof');
   lenminutes=ftell(f)/2/sr/60;%each sample is two bytes
   lenminutes=min(lenminutes, p.max_minutes2cluster);
   nchunks=ceil(lenminutes/handles.par.chunk);
   chunklenseconds=round(lenminutes/nchunks*60);
   begs=[0:nchunks-1]*chunklenseconds;
   ends=[[1:nchunks-1]*chunklenseconds lenminutes*60];
   
   spikes=[];
   index=[];
   fprintf('File splitted into %d chunks of %d seconds:\n',nchunks,chunklenseconds);
   for ich=1:nchunks,
      fprintf('%d ',ich);
      [spikes_, index_, thr]=Extract_spikes_ION(handles,begs(ich),ends(ich));
      fprintf('%d spikes detected, thr = %1.2f\n',size(spikes_,1),thr);
      spikes=[spikes; spikes_];
      index=[index; index_+begs(ich)*1000];
   end
   clear spikes_ index_
   if p.invert, spikes=-spikes;end
      
   if isempty(spikes), fprintf('No spikes detected in file %s\n',filename); continue; end
   nspk = size(spikes,1);
   if nspk<15, fprintf('Less than 15 spikes detected in file %s\n',filename); continue; end
   fprintf('%d spikes detected\n',nspk);
   
   handles.fname = sprintf('data_%s-%d',handles.bname,channel);   %filename for interaction with SPC

   min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
   handles.par.min_clus = min_clus;

   output_filename=sprintf('times_%s-%d.mat',handles.bname,channel+fmin);

   fprintf('Feature detection...\n');
   %CALCULATES INPUTS TO THE CLUSTERING ALGORITHM
   %get features
   [features,feature_names] = wave_features_ION(spikes,handles);
   save(output_filename,'features','feature_names');
 
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
   clear features feature_names

   
   [clu, tree] = run_cluster(handles);
   classtemp=cell(size(clu,1),handles.par.max_nrclasses);
   for te=1:size(clu,1),
      %       for j=1:handles.par.max_nrclasses, classtemp{te,j}=inds_to_cluster(find(clu(te,3:end)==j-1)); end;
      for j=1:handles.par.max_nrclasses, classtemp{te,j}=inds_to_cluster((clu(te,3:end)==j-1)); end;
   end
   fprintf('Adding SPC output information into results file...\n');
   %    save(output_filename,'tree','classtemp','-append');
   %    save(output_filename,'spikes','-append');
   save(output_filename,'tree','classtemp','spikes','-append');

   [temp] = find_temp(tree,handles);
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
   save(output_filename,'cluster_class','par','-append','-v6');

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
         eval(sprintf('print -djpeg fig2print_%s-%d',handles.bname,channel+fmin));
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
clear handles
