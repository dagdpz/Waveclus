function [inspk,feature_names,inputs] = wc_feature_selection(spikes,ts,handles)
%Calculates the spike features
scales = handles.WC.scales;
feature = handles.WC.features;
exclusioncrit = handles.WC.exclusioncrit; % thr; number
exclusionthr = handles.WC.exclusionthr;
maxinputs = handles.WC.maxinputs;
int_factor = handles.WC.int_factor;
w_pre = handles.WC.w_pre;
w_post = handles.WC.w_post;
nspk=size(spikes,1);
len = size(spikes,2)/int_factor;

%% this part is more than dubious...!!!!
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% w_pre=10;
% int_factor=1;
% shift = 2;
% w_post = 22;
offstart = round(w_pre/2) + shift;
offend  = round(w_post/2) + shift + w_pre;
ind1 = [fliplr(w_pre : -1: offstart+1),w_pre+1 : offend] * int_factor;

% CALCULATES FEATURES
ccall = [];
fnall = {};

%% Wavelet decomposition
if strfind(feature,'wav')
    cc=zeros(nspk,length(ind1));
    for i=1:nspk
        [c,l]=wavedec(spikes(i,ind1),scales,handles.WC.wavelet);
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
    norm_factor_WL=std(cc(:));
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
    S=S(:,1:numdim);
    norm_factor_PCA=std(S(:));
    S=S*norm_factor_WL/norm_factor_PCA;
    ccall = cat(2,ccall,S);
    fnall = cat(2,fnall,fn);
    clear fn S
end


%% Raw waveforms

if strfind(feature,'raw')
    spikesdown = spikes(:,ind1);
    %features
    fn = cell(1,length(ind1));
    for i=1:length(ind1)
        fn{i}=sprintf('Raw,%d',i);
    end
    norm_factor_raw=std(spikesdown(:));
    spikesdown=spikesdown*norm_factor_WL/norm_factor_raw;
    ccall = cat(2,ccall,spikesdown);
    fnall = cat(2,fnall,fn);
    clear fn inddummy spikesdown
end

%% time 
if strfind(feature,'time')
    norm_factor_t=1800000; %% normalized to 30 minutes
    draw=ts*norm_factor_WL/norm_factor_t;
    ccall = cat(2,ccall,draw);
    fnall = cat(2,fnall,{'Time'});
end

%% feature selection
warning('off','stats:lillietest:OutOfRangeP');
warning('off','stats:lillietest:OutOfRangePHigh');
warning('off','stats:lillietest:OutOfRangePLow');

stdcc = std(ccall) * 3;
meancc = mean(ccall);
thr_dist_min = meancc - stdcc;
thr_dist_max = meancc + stdcc;

sd=nan(1,size(ccall,2));
for i=1:size(ccall,2)                            % lilliefors test for coefficient selection
    aux = ccall(ccall(:,i) > thr_dist_min(i) & ccall(:,i) < thr_dist_max(i),i); % removing outliers (3x std)
    if length(aux) > 10;
        [~,~,lstat]=lillietest(aux,0.05);
        sd(i)=lstat;
    else
        sd(i)=0;
    end
end
if strfind(feature,'time')
sd(end)=max(sd)+1; %% to alway keep time??
end

[sortsd,ind]=sort(-sd); % orders indexes by lstat (KSstatistic) of lilliefors test, the higher ones first

% excluding features if KStatistics below certain threshold (mean +
% 1.96*std) of a bootstrapped random dummie
% if strcmp(exclusioncrit,'thr')
%     crit = nan(100,1);
%     for i = 1 : 100
%         dummy = randn(size(spikes,1),1);
%         dummy = dummy(dummy > -3 & dummy < 3);
%         [~,~,lstat]=lillietest(dummy,0.05);
%         crit(i) = lstat;
%     end
%     crit = mean(crit) + 1.96*std(crit);
%     
%     ind = ind(sortsd < -crit);
% end

warning('on','stats:lillietest:OutOfRangeP');
warning('on','stats:lillietest:OutOfRangePHigh');
warning('on','stats:lillietest:OutOfRangePLow');

ccall = ccall(:,ind);
fnall = fnall(ind);

clear ind

%% removing correlated features! we should probably be more strict here -> lower exclusionthr
if ~isempty(ccall)
    feat_corr = corr(ccall).^2;                    % feature correlation matrix
    feat_corr_true = feat_corr >= exclusionthr;           % logical matrix of correlated features
    feat_corr_true(diag(true(1,size(feat_corr_true,1)))) = 0; % removing autocorrelation
    ind = 1 : size(feat_corr_true,1);
    indrej = [];
    
    %% removing features one by one, starting with the one that correlates with most
    while sum(feat_corr_true(:))
        distdummy = sum(feat_corr_true);
        %thrdummy = exclusionthr;
        delind = find(distdummy == max(distdummy));
        if numel(delind) > 1 % several features with same max correlation
            simmat3 = feat_corr(delind,:);
            simmat3(simmat3<exclusionthr)=2;
            simmat3=round(simmat3*100)/100; % 0.01 stepsize maintained
            min_above_thr=min(simmat3,[],2);
            delind(find(min_above_thr==max(min_above_thr),1, 'last'));
        end
        feat_corr(:,delind) = [];
        feat_corr(delind,:) = [];
        feat_corr_true(:,delind) = [];
        feat_corr_true(delind,:) = [];
        indrej = cat(2,indrej,ind(delind));
        ind(delind)= [];
    end
    
    ccall = ccall(:,ind);
    fnall = fnall(ind);
end
if size(ccall,2) < maxinputs
    maxinputs = size(ccall,2);
end

inspk = ccall(:,1:maxinputs);
feature_names=fnall(1:maxinputs);
inputs = maxinputs;

fprintf('Selected features: ');
fprintf('%s ',feature_names{:});
fprintf('\n')
