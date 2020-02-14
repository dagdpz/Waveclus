
clear all
dbstop if error

spikecolors = ['k','b','r','g','c','m','y','b','r','g','c','m','y'];
scrsz = get(0,'ScreenSize');

clustofuse = xlsread('wrongsep');
clustofuse(:,1) = [];

% files = dir('times_datahighpass*.mat');
%files = dir('timesnew_datahighpass*.mat');
files = dir('dataspikes_*_redetect.mat');

filesB = {files.name};
for m= 1:length(filesB);
    tic
    filename = filesB{1,m};
    load(filename);
    
    cluster_classnew = cluster_class;
    dummy = clustofuse(str2double(filename(end-15:end-13)),:);
    
    if ~(sum(isnan(dummy)) == size(clustofuse,2))
        
        dummy2 = ~isnan(cat(2,nan,dummy,nan));
        dummy3 = diff(dummy2);
        indstart = find(dummy3 == 1);
        indend = find(dummy3 == -1)-1;
        
        
        for i = 1 : length(indstart)
            indtofuse = dummy(indstart(i):indend(i));
            for j = 1 : length(indtofuse)-1
                cluster_classnew(cluster_classnew(:,1) == indtofuse(j+1),1) = indtofuse(1);
            end
            clear indtofuse
        end
        
        clear dummy2 dummy3 indstart indend
        
        clusterind = unique(cluster_classnew(:,1));
        clusterind(clusterind == 0) = [];
        
        for i = 1 : length(clusterind)
            cluster_classnew(cluster_classnew(:,1) == clusterind(i),1) = i;
        end
        
%         diffcluster = setdiff(1:length(clusterind),clusterind);
%         if ~isempty(diffcluster)
%             for i = 1 : length(clusterind) - (diffcluster(1)-1)
%                 cluster_classnew(cluster_classnew(:,1) == ,1) = 
        clear clusterind        
    end
    
    clear dummy
    
    spikeall = [];
    ncl = max(cluster_classnew(:,1));
    
    for j = 1 : ncl
        spikeall = cat(1,spikeall,mean(spikes(cluster_classnew(:,1) == j,:)));
    end
    thr = 3*median(abs(cat(1,spikes(:,1),spikes(:,end))))/0.6745;
    ybounds = [min(spikeall(:))-thr max(spikeall(:))+thr];
    
    
    dummy = max(cluster_classnew(:,1));
    if dummy > 5
        cluster_classA = cluster_classnew(cluster_classnew(:,1) <= 5,:);
        spikesA = spikes(cluster_classnew(:,1) <= 5,:);
        spikeimg2(spikesA,cluster_classA,ybounds);
        set(gcf, 'PaperPositionMode', 'auto');
        print('-djpeg', ['datahighpassnew' filename(end-15:end-13) 'A'])
        close
        clear cluster_classA spikesA
        
        cluster_classB = cluster_classnew(cluster_classnew(:,1) > 5,:);
        spikesB = spikes(cluster_classnew(:,1) > 5,:);
        cluster_classB(:,1) = cluster_classB(:,1) - 5;
        spikeimg2(spikesB,cluster_classB,ybounds);
        set(gcf, 'PaperPositionMode', 'auto');
        print('-djpeg', ['datahighpassnew' filename(end-15:end-13) 'B'])
        close
        clear cluster_classB spikesB
        
    else
        spikeimg2(spikes,cluster_classnew,ybounds);
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpng', ['datahighpassnew' filename(end-15:end-13)])
        close
        
    end
    
    clear spikeall ncl thr ybounds
    
    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
    for i = 1 : max(cluster_classnew(:,1))
        
        spiketimes = zeros(round(cluster_classnew(end,2))*30,1);
        spiketimes(round(cluster_classnew(cluster_classnew(:,1) == i,2)*30)) = 1;
        spiketimes = spiketimes(1 : floor(length(spiketimes)/300000)*300000);
        spiketimes = sum(reshape(spiketimes,[300000 length(spiketimes)/300000]));
    
        subplot(max(cluster_classnew(:,1)),1,i);plot(spiketimes,spikecolors(i+1),'linewidth',2)
        clear spiketimes
    end
    print('-dpng', ['datahighpassnewhist' filename(end-15:end-13)])
    close
    
    save(filename,'cluster_classnew','-append');
    
    clear spikes cluster_class cluster_classnew
    toc
    
end