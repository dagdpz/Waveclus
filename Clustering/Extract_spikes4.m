function Extract_spikes4(handles)

sr=handles.par.sr;
w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
ref=handles.par.ref;
int_factor=handles.par.int_factor;

%files = dir('datafilt2_ch*.mat');
files = dir('datafilt_ch*.mat');
files = {files.name};

for k= 1:length(files)
    tic
    load(files{k},'data')
    
    data = double(data);
    %     if ~strcmp(handles.dtypewrite,'single')
    %         data=single(data);
    %     end
    
    if ~(handles.par.transform_factor == 1)
        data = data * handles.par.transform_factor;                       %%in microvolts
    end
    
    if size(data,2) == 1
        data = data';
    end
    
    numsamples = length(data);
    
    thr = handles.par.stdmin*median(abs(data))/0.6745;
    thrmax = handles.par.stdmax*median(abs(data))/0.6745;
    
    % LOCATE SPIKE TIMES
    switch handles.threshold
        case 'pos', ups = find(data(w_pre+2:end-w_post-2) > thr) +w_pre+1;%indices of values above threshold
        case 'neg', ups = find(data(w_pre+2:end-w_post-2) < -thr) +w_pre+1;%indices of values below threshold
        case 'both', ups = find(abs(data(w_pre+2:end-w_post-2)) > thr) +w_pre+1;
    end
    if isempty(ups), %no spikes detected
        spikes=[]; index=[];
        %    save(spikes_fname,'spikes','index');
    else
        %% this is the part to reduce to during the task only
        is_during_task=false(size(ups));
        for t=1:size(handles.task_times,1)
        is_during_task(ups>=handles.task_times(t,1)*handles.par.sr & ups<=handles.task_times(t,2)*handles.par.sr)=true;
        end
        ups=ups(is_during_task);
        %find beginnings and endings of intervals with values above threshold
        begs=ups([true diff(ups)>1]);%indices in data array
        ends=ups([diff(ups)>1 true]);
        clear ups is_during_task 
        
        max_index=numel(data)-w_post-2;
        min_index=w_pre+2;
        begs(begs<min_index)=min_index;
        begs(begs>max_index)=max_index;
        ends(ends<min_index)=min_index;
        ends(ends>max_index)=max_index;
        
        nspk=length(begs);
        index=nan(nspk,1);
        direction=nan(nspk,1);
        np=w_post+w_pre;
        spikes=nan(nspk,np+4);
        m=nan(nspk,1);
        for i=1:nspk,
            [m(i),ind] = max(abs(data(begs(i):ends(i))));
            index(i) = begs(i) + ind -1;
            spikes(i,:)= data(index(i)-w_pre+1-2:index(i)+w_post+2);
            %% modification here LS 20171010
            direction(i)=sign(data(index(i)));
        end        
        % Eliminate artifacts
        index(m > thrmax)=[];
        spikes(m > thrmax,:)=[];
        direction(m > thrmax,:)=[];
        m(m > thrmax,:)=[];
        
        clear data begs ends %Tesla 20160608: Out of memory
        
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
            %% LS 20171010 this is the modification for removing positive spikes in case there was both
            directionality=direction(fff);
            neg_prefindex=any(directionality<0,2) & any(directionality>0,2);
            [fff_idx(neg_prefindex),~]=find(directionality(neg_prefindex,:)'>0);
            
            fff=[fff(fff_idx==1,1);fff(fff_idx==2,2)];
            direction(fff)=[];
            index(fff)=[];
            m(fff)=[];
            spikes(fff,:)=[];
        end
        
        nspk=numel(index);
        
        %INTERPOLATION
        if handles.par.interpolation=='y',
            %Does interpolation
            ints=1/int_factor:1/int_factor:np+4;
            intspikes = spline(1:np+4,spikes,ints);
            spikes1=zeros(nspk,np*int_factor);
            
            switch handles.threshold
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
        
        index=index/sr*1000;
        par=handles.par;
        cluster_class=[zeros(size(index)) index];
        %         chanceevent = numsamples * chancelevel/2;
        
        switch handles.threshold
            case 'pos'
                save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr','par','cluster_class')
           case 'neg'
                save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr','par','cluster_class')
            case 'both'
                spikesign = sign(spikes(:,w_pre*int_factor));
                spikesall = spikes;
                indexall = index;
                spikes = spikesall(spikesign == -1,:);
                index = indexall(spikesign == -1);
                save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_negthr'],'spikes','index','thr','par','cluster_class')
                spikes = spikesall(spikesign == 1,:);
                index = indexall(spikesign == 1);
                save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_posthr'],'spikes','index','thr','par','cluster_class')
         end
        
        
    end
    clearvars -except handles sr ref w_pre w_post int_factor files k chancelevel
    toc
end


