function wc_extract_spikes_cat_MU_SU(handles)

sr=handles.WC.sr;
w_pre=handles.WC.w_pre;
w_post=handles.WC.w_post;
ref=handles.WC.ref;
int_factor=handles.WC.int_factor;
remove_ini=handles.WC.remove_ini;

%% entirely different loop here... only loop to concatinate data, outside we loop through channels (&blocks)
dat=[];
for b=1:numel(handles.current_blocks)
    load([handles.WC_concatenation_folder(1:end-3) 'WC_Block-' num2str(handles.current_blocks(b)) filesep 'datafilt_ch'  sprintf('%03d',handles.current_channel) '.mat'],'data');
    dat=[dat data];
    clear data
end



tic
%% this might be a source for problems ... apparently not?
%     dat = double(dat);
%     if ~strcmp(handles.dtypewrite,'single')
%         dat=single(dat);
%     end

if ~(handles.WC.transform_factor == 1)
    dat = dat * handles.WC.transform_factor;                       %%in microvolts
end

if size(dat,2) == 1
    dat = dat';
end

%numsamples = length(dat);

thr = handles.WC.stdmin*median(abs(dat))/0.6745;
thrmax = handles.WC.stdmax*median(abs(dat))/0.6745;

% LOCATE SPIKE TIMES
switch handles.WC.threshold
    case 'pos', ups = find(dat(w_pre+2:end-w_post-2) > thr) +w_pre+1;%indices of values above threshold
    case 'neg', ups = find(dat(w_pre+2:end-w_post-2) < -thr) +w_pre+1;%indices of values below threshold
    case 'both', ups = find(abs(dat(w_pre+2:end-w_post-2)) > thr) +w_pre+1;
end
    %% this is the part to reduce to during the task only
    if remove_ini
        is_during_task=false(size(ups));
        for t=1:size(handles.task_times,1)
            is_during_task(ups>handles.task_times(t,1)*handles.WC.sr & ups<handles.task_times(t,2)*handles.WC.sr)=true; %% changed from <= and >= to </>
        end
        ups=ups(is_during_task);
    end
if isempty(ups), %no spikes detected
    spikes=[]; index=[]; spikesign=[];
    %    save(spikes_fname,'spikes','index');
else
    %find beginnings and endings of intervals with values above threshold
    begs=ups([true diff(ups)>1]);%indices in dat array
    ends=ups([diff(ups)>1 true]);
    clear ups is_during_task
    
    max_index=numel(dat)-w_post-2;
    min_index=w_pre+2;
    begs(begs<min_index)=min_index;
    begs(begs>max_index)=max_index;
    ends(ends<min_index)=min_index;
    ends(ends>max_index)=max_index;
    
    nspk=length(begs);
    index=nan(nspk,1);
    spikesign=nan(nspk,1);
    np=w_post+w_pre;
    spikes=nan(nspk,np+4);
    m=nan(nspk,1);
    for i=1:nspk,
        [m(i),ind] = max(abs(dat(begs(i):ends(i))));
        index(i) = begs(i) + ind -1;
        spikes(i,:)= dat(index(i)-w_pre+1-2:index(i)+w_post+2);
        %% modification here LS 20171010
        spikesign(i)=sign(dat(index(i)));
    end
    % Eliminate artifacts
    index(m > thrmax)=[];
    spikes(m > thrmax,:)=[];
    spikesign(m > thrmax,:)=[];
    m(m > thrmax,:)=[];
    
    clear dat begs ends %Tesla 20160608: Out of memory
    
    % Eliminate spikes with too short refractory period.
    % Keep larger amplitude BUT prefer negative part LS 20171010
    while true
        fff=find(diff(index)<=ref*sr);
        if isempty(fff)
            break
        end
        fff=[fff fff+1];
        if numel(fff)==2
            [~, fff_idx]=min(m(fff),[],1);
        else
            [~, fff_idx]=min(m(fff),[],2);
        end
%         %% LS 20171010 this is the modification for removing positive spikes in case there was both
%         directionality=spikesign(fff);
%         neg_prefindex=any(directionality<0,2) & any(directionality>0,2);
%         [fff_idx(neg_prefindex),~]=find(directionality(neg_prefindex,:)'>0);
        
        fff=[fff(fff_idx==1,1);fff(fff_idx==2,2)];
        spikesign(fff)=[];
        index(fff)=[];
        m(fff)=[];
        spikes(fff,:)=[];
    end
    
    nspk=numel(index);
    
    %INTERPOLATION
    if handles.WC.interpolation=='y',
        %Does interpolation
        ints=1/int_factor:1/int_factor:np+4;
        intspikes = spline(1:np+4,spikes,ints);
        spikes1=zeros(nspk,np*int_factor);
        
        switch handles.WC.threshold
            case 'pos'
                [~,iaux]=max(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
            case 'neg'
                [~,iaux]=min(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
            case 'both'
                [~,iaux]=max(abs(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor)),[],2);
        end
        iaux = iaux + (w_pre+1)*int_factor -1;
        for i=1:nspk
            spikes1(i,w_pre*int_factor:-1:1) = intspikes(i,iaux(i):-1:iaux(i)-(w_pre*int_factor-1));
            spikes1(i,w_pre*int_factor+1:np*int_factor) = intspikes(i,iaux(i)+1:1:iaux(i)+w_post*int_factor);
        end
        spikes=spikes1;
        clear spikes1
        clear intspikes
    else
        spikes(:,[1 2 end-1 end])=[];       %eliminates borders that were introduced for interpolation
    end
end

index=index/sr*1000;
par=handles.WC;
cluster_class=[zeros(size(index)) index];
fnamepart=[handles.WC_concatenation_folder 'dataspikes_ch' sprintf('%03d',handles.current_channel) '_' num2str(handles.current_channel_file) '_' handles.current_threshold_step];
switch handles.WC.threshold
    case 'pos'
        save([fnamepart '_pos'],'spikes','index','thr','par','cluster_class')
    case 'neg'
        save([fnamepart '_neg'],'spikes','index','thr','par','cluster_class')
    case 'both'
        %spikesign = sign(spikes(:,w_pre*int_factor));
        spikesall = spikes;
        indexall = index;
        cluster_class_all=cluster_class;
        spikes = spikesall(spikesign == -1,:);
        cluster_class = cluster_class_all(spikesign == -1,:);
        index = indexall(spikesign == -1);
        save([fnamepart '_neg'],'spikes','index','thr','par','cluster_class')
        spikes = spikesall(spikesign == 1,:);
        cluster_class = cluster_class_all(spikesign == 1,:);
        index = indexall(spikesign == 1);
        save([fnamepart '_pos'],'spikes','index','thr','par','cluster_class')
end

clearvars -except handles sr ref w_pre w_post int_factor files k chancelevel
    disp([' Extracting spikes ' handles.WC_concatenation_folder(1:end-3) ' ch '  sprintf('%03d',handles.current_channel) ' took ' num2str(round(toc*10)/10) ' seconds']);       


