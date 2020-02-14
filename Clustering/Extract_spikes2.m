function Extract_spikes2(handles)

sr=handles.par.sr;
w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
int_factor=handles.par.int_factor;

files = textread('Files.txt','%s');

for k= 1:length(files)
    tic
    file_to_cluster = files(k)
    eval(['load ' char(file_to_cluster) ';']);
    
    x=double(data') * handles.par.transform_factor;  %%in microvolts
    clear data dataold
    
    thr = handles.par.stdmin*median(abs(x))/0.6745;
    thrmax = handles.par.stdmax*median(abs(x))/0.6745;
    
    % LOCATE SPIKE TIMES
    switch handles.threshold
        case 'pos', ups = find(x(w_pre+2:end-w_post-2) > thr) +w_pre+1;%indices of values above threshold
        case 'neg', ups = find(x(w_pre+2:end-w_post-2) < -thr) +w_pre+1;%indices of values below threshold
        case 'both', ups = find(abs(x(w_pre+2:end-w_post-2)) > thr) +w_pre+1;   
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
        %find two events which are closer than 0.2ms, remove smallest from two
        %closest
        bad=[];
        for i=1:length(begs)-1,
            if begs(i+1)-ends(i)<=1.4e-3*sr, bad=[bad i]; end
        end
        bad2=[];
        for i=1:length(begs)-2,
            if begs(i+2)-ends(i)<=1.4e-3*sr, bad2=[bad2 i]; end
        end

        nspk=length(begs);
        index=NaN*ones(nspk,1);
        np=w_post+w_pre;
        spikes=NaN*ones(nspk,np+4);
        m=NaN*ones(nspk,1);
        for i=1:nspk,
            [m(i),ind] = max(abs(x(begs(i):ends(i))));
            index(i) = begs(i) + ind -1;
            spikes(i,:)= x (index(i)-w_pre+1-2:index(i)+w_post+2);
        end
        clear x
        bad1=[];
        for i=bad,
            if m(i)>m(i+1), bad1=[bad1 i+1];
            else bad1=[bad1 i]; end
        end
        if ~isempty(bad2)  
        for i=bad2,
            if m(i)>m(i+2), bad1=[bad1 i+2];
            else bad1=[bad1 i]; end
        end
        end
         bad1 = unique(sort(bad1));
        ind=setdiff(1:nspk,bad1);
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
                    [maxi,iaux]=max(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
                case 'neg'
                    [maxi,iaux]=min(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
                case 'both'
                    [maxi,iaux]=max(abs(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor)),[],2);
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
            spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
            spikes(:,1:2)=[];
        end
        index=index/sr*1000;
        % save(spikes_fname,'spikes','index');
        eval(['save ' char(file_to_cluster) '_spikes.mat spikes index thr']);
    end
    clearvars -except handles sr w_pre w_post int_factor files k
    toc
end


