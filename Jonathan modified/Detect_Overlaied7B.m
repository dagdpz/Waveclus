clear all
dbstop if error

spikecolors = ['c';'b';'r';'g';'y';'m'];
scrsz = get(0,'ScreenSize');

handles.classify_space = 'spikeshapes';
% handles.classify_method = 'diaglinear';
handles.classify_method = 'linear';
handles.threshold ='both';

handles.sdnum = 5;
repetitions = 2;

chunks = 3000000;

% tempmeth = 'Tempmatch';
tempmeth = 'Classify';
int_factor = 2;

% % channels = textread('Files.txt','%s');
% channels = dir('datahighpass*.mat');
% channels = {channels.name};

spikefilenames = 'dataspikes_ch';

files = dir([spikefilenames '*.mat']);
filesB = {files.name};

chan2redetect = cell(length(filesB),1);
for i = 1 : length(filesB)
    chan2redetect{i} = filesB{i}(1:length(spikefilenames)+3);
end
chan2redetect = unique(chan2redetect);

    
% cellfun(@mat2double, filesB
% save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_negthr'],'spikes','index','thr')
% str2
% [] = unique


for m = 1:length(chan2redetect)
    tic
    filename = chan2redetect{m};
    disp(filename)
    spikesAll = [];
    cluster_classAll = [];
    numofclus = 0;
    if exist([filename '_negthr.mat'],'file')
        load([filename '_negthr.mat'],'spikes','cluster_class','par','thr');
        spikesAll = cat(1,spikesAll,spikes);
        cluster_classAll = cat(1,cluster_classAll,cluster_class);
        numofclus = max(cluster_class(:,1));
        clear spikes cluster_class
    end
    if exist([filename '_posthr.mat'],'file')
        load([filename '_posthr.mat'],'spikes','cluster_class','par','thr');
        cluster_class(cluster_class(:,1) > 0,1) = cluster_class(cluster_class(:,1) > 0,1) + numofclus;
        spikesAll = cat(1,spikesAll,spikes);
        cluster_classAll = cat(1,cluster_classAll,cluster_class);
        clear spikes cluster_class
    end
    [~, ind] = sort(cluster_classAll(:,2));
    spikes = spikesAll(ind,:);
    cluster_class = cluster_classAll(ind,:);
    clear ind spikesAll cluster_classAll
    
    handles.par = par;
%     if isempty(handles.par.thr)
%         load([filename(7:end-4) '_spikes.mat'],'thr')
%         handles.par.thr = thr;
%         clear thr
%     end 
    handles.par.thr = thr;
    clear thr
    delartrange = -0.2e-3*handles.par.sr: 0.2e-3*handles.par.sr;
    shift = 2;
    shift = shift*handles.par.int_factor;
    spikesold = spikes;
    clear spikes par features
    handles.ncl = max(cluster_class(:,1));
    handles.nspk = size(spikesold,1);
    spikeall = [];
    maxdist = [];
    pointdist = [];
    for j = 1 : handles.ncl
        spikeall = cat(1,spikeall,mean(spikesold(cluster_class(:,1) == j,:)));
        %         maxdist = cat(1,maxdist,sqrt(sum(var(spikes(cluster_class(:,1) == j,:),1))));
        %         pointdist   = cat(1,pointdist,sqrt(var(spikes(cluster_class(:,1) == j,:),1)));
        handles.classind{j} = find(cluster_class(:,1) == j)';
    end
    
    spikeall2 = [];
    spikeedge = abs(spikeall(:,1)) > handles.par.thr/handles.par.stdmin*2;
    if sum(spikeedge) 
        load(filename(7:end));
        dummy = 1:size(spikeall,1);
        spikeall2 = nan(size(spikeall,1),size(spikeall,2)+handles.par.w_pre*2);
        for k = dummy(spikeedge)
            index = round(cluster_class(cluster_class(:,1) == k,2)*(handles.par.sr/1000));
            nspk = size(index,1);
            
            spikes = data(repmat(-(handles.par.w_pre*2)+1-2:handles.par.w_post+2,[nspk 1])+repmat(index,[1 handles.par.w_pre*2+handles.par.w_post+4]));
            
            if handles.par.interpolation=='y',
                %Does interpolation
                ints=1/int_factor:1/int_factor:(handles.par.w_post+2*handles.par.w_pre)+4;
                intspikes = spline(1:(handles.par.w_post+2*handles.par.w_pre)+4,double(spikes),ints);
                clear spikes
                spikes=zeros(nspk,(handles.par.w_post+2*handles.par.w_pre)*int_factor);
                
                switch handles.threshold
                    case 'pos'
                        [maxi,iaux]=max(intspikes(:,(handles.par.w_pre*2+1)*int_factor:(handles.par.w_pre*2+3)*int_factor),[],2);
                    case 'neg'
                        [maxi,iaux]=min(intspikes(:,(handles.par.w_pre*2+1)*int_factor:(handles.par.w_pre*2+3)*int_factor),[],2);
                    case 'both'
                        [maxi,iaux]=max(abs(intspikes(:,(handles.par.w_pre*2+1)*int_factor:(handles.par.w_pre*2+3)*int_factor)),[],2);
                end
                iaux = iaux + (handles.par.w_pre*2+1)*int_factor -1;
                for i=1:nspk
                    spikes(i,handles.par.w_pre*2*int_factor:-1:1) = intspikes(i,iaux(i):-1:iaux(i)-(handles.par.w_pre*2*int_factor-1));
                    spikes(i,handles.par.w_pre*2*int_factor+1:(handles.par.w_post+2*handles.par.w_pre)*int_factor) = intspikes(i,iaux(i)+1:1:iaux(i)+handles.par.w_post*int_factor);
                end
                
                
                spikeall2(k,:) = mean(spikes);
                clear intspikes iaux spikes
            end
        end
        clear data dummy index nspk
    end
    
    clear cluster_class
    
    spikediffal = nan(size(spikeall,1),size(spikeall,2)/2,3);
    
    spikediffal(:,:,1) = spikeall(:,1:2:end-1);
    spikediffal(:,:,2) = spikeall(:,2:2:end);
    spikediffal(:,:,3) = cat(2,spikeall(:,3:2:end),spikeall(:,end));
    
    if ~isempty(spikeall2)
    spikediffal2 = nan(size(spikeall2,1),size(spikeall2,2)/2,3);
    
    spikediffal2(:,:,1) = spikeall2(:,1:2:end-1);
    spikediffal2(:,:,2) = spikeall2(:,2:2:end);
    spikediffal2(:,:,3) = cat(2,spikeall2(:,3:2:end),spikeall2(:,end));
    end
    
    [~,indexspikemax]=sort(-abs(spikeall(:,handles.par.w_pre*int_factor)));
    
    classes=zeros(1,handles.nspk);
    for i=1:handles.ncl
        classes(handles.classind{i})=i;
    end
    
    handles.classind = [];
    
    group=classes([classes~=0]);
    spikesold(classes == 0,:) = [];
    clear classes
    
    ints1 = sort(handles.par.w_pre*handles.par.int_factor: -handles.par.int_factor: handles.par.w_pre*handles.par.int_factor-handles.par.w_pre*handles.par.int_factor/2+1+shift);
    ints2 = handles.par.w_pre*handles.par.int_factor + handles.par.int_factor : handles.par.int_factor: handles.par.w_pre*handles.par.int_factor+handles.par.w_post*handles.par.int_factor/2+shift;
    handles.ints = cat(2,ints1,ints2);
    clear ints1 ints2
    
%     channel=channels{m};
%     load(channel)
    load(['datafilt2' filename(end-5:end)])
    clear dataold
    data = double(data);
    cluster_classpre = [];
    spikespre = [];
    thrall = [];
    
    for repfac = 1 : repetitions
        
        for unit = 1 : size(indexspikemax,1)
            
            if repfac == 1;
                throld = [];
            else
                throld = thrall(unit);
            end
            
            training = cat(2,spikesold(:,handles.ints));
            dataconv = conv(data,-fliplr(spikeall(indexspikemax(unit),int_factor:int_factor:end)));
            dataconv(end-size(spikeall,2)+2 : end) = [];
            
            [spikes,index,thr,subshift] = Extract_spikes_Overlaiedopfi6(dataconv,data,repfac,throld,handles.par);
            if repfac == 1
                thrall = cat(1,thrall,thr);
            end
            
            clear dataconv thr throld
            sample=cat(2,spikes(:,handles.ints));
            %             clear spikes
            
            switch tempmeth
                case 'Tempmatch'
                    class_out = zeros(1,size(cluster_class,1));
                    for k=1:nspk,
                        class_out(k) = nearest_neighbor(spikes(k,:),newspikes,sdnum*maxdist);
                    end
                    
                case 'Classify'
                    handles=classifyrestoverlayed3(handles,training,sample,group);
                    sample = [];
            end
            
            Unittimes = index(handles.c == indexspikemax(unit));
            Unitshift = subshift(handles.c == indexspikemax(unit));
            
            spikes2 = spikes(handles.c == indexspikemax(unit),:);
            if unit * repfac > 1
                [~,indreal] = setdiff(round(Unittimes*(handles.par.sr/1000)),reshape(repmat(delartrange,[size(cluster_classpre,1) 1]) + repmat(round(sort(cluster_classpre(:,2))*(handles.par.sr/1000)),[1 length(delartrange)]),[size(cluster_classpre,1)*length(delartrange) 1]));
                Unittimes = Unittimes(indreal);
                Unitshift = Unitshift(indreal);
                spikes2 = spikes2(indreal,:);
                clear indreal
            end
            
            if spikeedge(indexspikemax(unit))
                for j = 1 : length(Unittimes)
                    if Unitshift(j) == 2
                        data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal2(indexspikemax(unit),:,3))';
                    elseif Unitshift(j) == 3
                        data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal2(indexspikemax(unit),:,2))';
                    elseif Unitshift(j) == 4
                        data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-40+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal2(indexspikemax(unit),:,1))';
                    end
                end
                
            else
                
                for j = 1 : length(Unittimes)
                    if Unitshift(j) == 2
                        data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal(indexspikemax(unit),:,3)/handles.par.transform_factor);
                    elseif Unitshift(j) == 3
                        data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal(indexspikemax(unit),:,2)/handles.par.transform_factor);
                    elseif Unitshift(j) == 4
                        data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44) = data(round(Unittimes(j) * (handles.par.sr/1000))-20+1:round(Unittimes(j) * (handles.par.sr/1000))+44)-(spikediffal(indexspikemax(unit),:,1)/handles.par.transform_factor);
                    end
                    
                end
            end
            
            if ~isempty(Unittimes)
                cluster_classpre = cat(1,cluster_classpre,cat(2,repmat(indexspikemax(unit),[length(Unittimes) 1]),Unittimes));
                
                spikespre = cat(1,spikespre,spikes2);
            end
            clear spikes2 Unittimes
            
            if unit*repfac == size(indexspikemax,1)*repetitions
                Unittimes = index(handles.c == 0);
                spikes2 = spikes(handles.c == 0,:);
                [~,indreal] = setdiff(round(Unittimes*(handles.par.sr/1000)),reshape(repmat(delartrange,[size(cluster_classpre,1) 1]) + repmat(round(sort(cluster_classpre(:,2))*(handles.par.sr/1000)),[1 length(delartrange)]),[size(cluster_classpre,1)*length(delartrange) 1]));
                Unittimes = Unittimes(indreal);
%                 Unitshift = Unitshift(indreal);
                spikes2 = spikes2(indreal,:);
                clear indreal
                if ~isempty(Unittimes)
                    cluster_classpre = cat(1,cluster_classpre,cat(2,zeros(length(Unittimes),1),Unittimes));
                    spikespre = cat(1,spikespre,spikes2);
                end
            end
            
            clear spikes training spikesadd spikescomp index subshift Unittimes Unitshift spikes2
            handles.c = [];
        end
    end
    
    thr = 3*median(abs(data))/0.6745;
    ybounds = [min(spikeall(:))-thr max(spikeall(:))+thr];
    
    clear indexspikemax data spikediffal spikesaddold spikescompold spikesold group spikeall thrall spikeall2
%     group = [];
    handles.par = [];
    
    
    [~,IX] = sort(cluster_classpre(:,2));
    spikes = spikespre(IX,:);
    cluster_class = cluster_classpre(IX,:);
    
    clear spikespre cluster_classpre IX

    dummy = max(cluster_class(:,1));
    if dummy > 5
    cluster_classA = cluster_class(cluster_class(:,1) <= 5,:);
    spikesA = spikes(cluster_class(:,1) <= 5,:);  
    spikeimg2(spikesA,cluster_classA,ybounds);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-djpeg', ['datahighpass' filename(14:end) 'A'])
    close
    clear cluster_classA spikesA 
    
    cluster_classB = cluster_class(cluster_class(:,1) > 5,:);
    spikesB = spikes(cluster_class(:,1) > 5,:); 
    cluster_classB(:,1) = cluster_classB(:,1) - 5;
    spikeimg2(spikesB,cluster_classB,ybounds);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-djpeg', ['datahighpass' filename(14:end) 'B'])
    close
    clear cluster_classB spikesB
    
    else
    spikeimg2(spikes,cluster_class,ybounds);  
    set(gcf, 'PaperPositionMode', 'auto');
    print('-djpeg', ['datahighpass' filename(14:end)])
    close   
        
    end
    
%     save([filename(1:5) '_datahighpass' filename(14:end)],'cluster_class','spikes');
    save([filename '_redetect'],'cluster_class','spikes');
    
    clear spikes cluster_class 
    toc
end
    