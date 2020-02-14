function Extract_spikes3(handles)

sr=handles.par.sr;
w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
int_factor=handles.par.int_factor;

% switch handles.par.stdmin
%     case 2
%         chancelevel = 0.0455;
%     case 3
%         chancelevel = 0.0027;
%     case 4
%         chancelevel = 0.000063;
%     case 5
%         chancelevel = 0.00000057;
% end
if handles.arraynoisecancelation
 files = dir('datafilt2_ch*.mat');
else
 files = dir('datafilt_ch*.mat');
end
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
        
        %find beginnings and endings of intervals with values above threshold
        begs=[0 find(diff(ups)>1)]+1;%index in index array ups
        ends=[begs(2:end)-1 length(ups)];
        begs=ups(begs);%indices in data array
        ends=ups(ends);
        %find two events which are closer than 1.4ms, remove smallest from two
        %closest
        bad=[];
        for i=1:length(begs)-1,
            if begs(i+1)-ends(i)<=1.4e-3*sr 
                bad=[bad i]; 
            end
        end
        bad2=[];
        for i=1:length(begs)-2,
            if begs(i+2)-ends(i)<=1.4e-3*sr 
                bad2=[bad2 i]; 
            end
        end

        nspk=length(begs);
        index=nan(nspk,1);
        np=w_post+w_pre;
        spikes=nan(nspk,np+4);
        m=nan(nspk,1);
        
        for i=1:nspk,
            [m(i),ind] = max(abs(data(begs(i):ends(i))));
            index(i) = begs(i) + ind -1;
            spikes(i,:)= data(index(i)-w_pre+1-2:index(i)+w_post+2);
        end
        clear data
        badind=[];
        if ~isempty(bad)
            for i=bad,
                if m(i)>m(i+1)
                    badind=[badind i+1];
                else
                    badind=[badind i];
                end
            end
        end
        if ~isempty(bad2)
            for i=bad2,
                if m(i)>m(i+2)
                    badind=[badind i+2];
                else
                    badind=[badind i];
                end
            end
        end
        badind = unique(sort(badind));
        ind=setdiff(1:nspk,badind);
        ind = ind(m(ind) <= thrmax);
        spikes=spikes(ind,:);
        index=index(ind);
        nspk=length(ind);
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
        
%         chanceevent = numsamples * chancelevel/2;    
        
        switch handles.threshold
            case 'pos'
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr')
%                 end
            case 'neg'
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr')
%                 end
            case 'both'
                spikesign = sign(spikes(:,w_pre*int_factor));
                spikesall = spikes;
                indexall = index;
                spikes = spikesall(spikesign == -1,:);
                index = indexall(spikesign == -1);
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_negthr'],'spikes','index','thr')
%                 end
                spikes = spikesall(spikesign == 1,:);
                index = indexall(spikesign == 1);
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_posthr'],'spikes','index','thr')
%                 end
        end
        
       
    end
    clearvars -except handles sr w_pre w_post int_factor files k chancelevel
    toc
end


