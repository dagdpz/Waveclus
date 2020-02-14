function [inspk,feature_names,inputs] = wave_features5(spikes,handles)
%Calculates the spike features

scales = handles.par.scales;
feature = handles.par.features;
selectionmeth = handles.par.selectionmeth; %eu: eucleadiean; corr: correlation; lilli: lillifors  
exclusioncrit = handles.par.exclusioncrit; % thr; number 
exclusionthr = handles.par.exclusionthr;
inputs = handles.par.inputs;
int_factor = handles.par.int_factor;
w_pre = handles.par.w_pre;
w_post = handles.par.w_post;
% sr = handles.par.sr;

nspk=size(spikes,1);
len = size(spikes,2)/int_factor;


shift = 2;
offstart = round(w_pre/2) + shift;
offend  = round(w_post/2) + shift + w_pre;

ind1 = [fliplr(w_pre : -1: offstart+1),w_pre+1 : offend] * int_factor;

% downfac = 2;
% ind2 = [fliplr(w_pre : -downfac : offstart+1),w_pre+1 : downfac : offend] * int_factor;

% CALCULATES FEATURES

ccall = [];
fnall = {};

%% Wavelet decomposition
if strfind(feature,'wav')
    
    cc=zeros(nspk,length(ind1));
    for i=1:nspk                                
        [c,l]=wavedec(spikes(i,ind1),scales,handles.par.wavelet);
        cc(i,:)=c; 
    end
    
    %features
    fn=cell(1,len/2);
    j=1;
    k=scales+1;
    ll=cumsum(l);
    for i=ll(1:end-1),
        j1=1;
        while j<=i,
            if i==ll(1), fn{j}=sprintf('A%d,%d',k-1,j1);
            else fn{j}=sprintf('D%d,%d',k,j1); end
            j=j+1; j1=j1+1;
        end
        k=k-1;
    end
    
    ccall = cat(2,ccall,cc);
    fnall = cat(2,fnall,fn);
    clear fn cc
end

%% PCA

if strfind(feature,'pca')
    numdim = 10;
    
    [~,S] = princomp(spikes(:,ind1));
    
    %features
    fn = cell(1,numdim);
    for i=1:numdim
        fn{i}=sprintf('PCA,%d',i);
    end
    
    ccall = cat(2,ccall,S(:,1:numdim));
    fnall = cat(2,fnall,fn);
    clear fn S
end


%% Raw waveforms

if strfind(feature,'raw')
    inddummy = ind1;
    inddummy(ind1 == w_pre * int_factor) = [];
    spikesdown = spikes(:,inddummy);
    
    %features
    fn = cell(1,length(ind1));
    for i=1:length(ind1)  
        fn{i}=sprintf('Raw,%d',i);
    end
    
    fn(ind1 == w_pre * int_factor) = [];
    
    ccall = cat(2,ccall,spikesdown);
    fnall = cat(2,fnall,fn);
    clear fn inddummy spikesdown
end

%% first order derivative of waveforms

if strfind(feature,'deriv')
       inddummy = ind1(1:round(length(ind1)/2)+1);
       draw = diff(spikes(:,inddummy),1,2);

       %features
    fn = cell(1,length(inddummy));
    for i=1:length(inddummy)
        fn{i}=sprintf('Deriv,%d',i);
    end
    
    ccall = cat(2,ccall,draw);
    fnall = cat(2,fnall,fn);
    clear fn inddummy draw
       
end


       
%% feature selection

switch selectionmeth
    case 'corr'
        
        simmat = corr(ccall).^2;
        simmat2 = simmat >= exclusionthr;
        simmat2(diag(true(1,size(simmat2,1)))) = 0;
        ind = 1 : size(simmat2,1);
        indrej = [];
        
        while sum(simmat2(:))
            distdummy = sum(simmat2);
            thrdummy = exclusionthr;
            delind = find(distdummy == max(distdummy));
            while length(delind) > 1
                simmat3 = simmat >= thrdummy;
                simmat3(diag(true(1,size(simmat3,1)))) = 0;
                distdummy = sum(simmat3);
                delind = delind(distdummy(delind) == max(distdummy(delind)));
                thrdummy = thrdummy+0.01;
                if ~sum(distdummy(delind))
                    delind = randsample(delind,1);
                end
            end
            simmat(:,delind) = [];
            simmat(delind,:) = [];
            simmat2(:,delind) = [];
            simmat2(delind,:) = [];
            indrej = cat(2,indrej,ind(delind));
            ind(delind)= [];
        end
        ind2 = cat(2,ind,fliplr(indrej));
        
    case 'eu'
        simmat = squareform(pdist(ccall','seuclidean')/sqrt(size(ccall,1)));
        simmat2 = simmat <= exclusionthr; %sqrt(2) + sqrt(1/size(ccall,1))*1.96; % correspond to 95% conf intervall for standart normal distributed noise difference
        simmat2(diag(true(1,size(simmat2,1)))) = 0;
        ind = 1 : size(simmat2,1);
        indrej = [];
        
        while sum(simmat2(:))
            distdummy = sum(simmat2);
            thrdummy = exclusionthr;
            delind = find(distdummy == max(distdummy));
            while length(delind) > 1
                simmat3 = simmat <= thrdummy;
                simmat3(diag(true(1,size(simmat3,1)))) = 0;
                distdummy = sum(simmat3);
                delind = delind(distdummy(delind) == max(distdummy(delind)));
                thrdummy = thrdummy*0.9;
                if ~sum(distdummy(delind))
                    delind = randsample(delind,1);
                end
            end
            simmat(:,delind) = [];
            simmat(delind,:) = [];
            simmat2(:,delind) = [];
            simmat2(delind,:) = [];
            indrej = cat(2,indrej,ind(delind));
            ind(delind)= [];
        end
        ind2 = cat(2,ind,fliplr(indrej));
        
    case 'lilli'
        
        warning('off','stats:lillietest:OutOfRangeP');
        warning('off','stats:lillietest:OutOfRangePHigh');
        warning('off','stats:lillietest:OutOfRangePLow');
        
        stdcc = std(ccall) * 3;
        meancc = mean(ccall);
        thr_dist_min = meancc - stdcc;
        thr_dist_max = meancc + stdcc;
        
        sd=nan(1,size(ccall,2));
        for i=1:size(ccall,2)                            % lilliefors test for coefficient selection
            aux = ccall(ccall(:,i) > thr_dist_min(i) & ccall(:,i) < thr_dist_max(i),i);
            if length(aux) > 10;
                [~,~,lstat]=lillietest(aux);
                sd(i)=lstat;
            else
                sd(i)=0;
            end
        end
        
        [sortsd,ind2]=sort(-sd);
        
        
        if strcmp(exclusioncrit,'thr')
            crit = nan(100,1);
            for i = 1 : 100
                dummy = randn(size(spikes,1),1);
                dummy = dummy(dummy > -3 & dummy < 3);
                [~,~,lstat]=lillietest(dummy);
                crit(i) = lstat;
            end
            crit = mean(crit) + 1.96*std(crit);
            
            ind = ind2(sortsd < -max(exclusionthr,crit));
          
        end
        
        warning('on','stats:lillietest:OutOfRangeP');
        warning('on','stats:lillietest:OutOfRangePHigh');
        warning('on','stats:lillietest:OutOfRangePLow');
                  
end

switch exclusioncrit
    case 'number'
        inspk = ccall(:,ind2(1:inputs));
        feature_names={fnall{ind2(1:inputs)}};
    case 'thr'    
       inspk = ccall(:,ind);
       feature_names={fnall{ind}};
       inputs = length(ind);
end

fprintf('Selected freatures: ');
fprintf('%s ',feature_names{:});
fprintf('\n')

% %CREATES INPUT MATRIX FOR SPC
% inspk=zeros(nspk,inputs);
% for i=1:nspk
%     for j=1:inputs
%         inspk(i,j)=cc(i,coeff(j));
%     end
% end

