function [inspk,feature_names,inputs] = wc_feature_selection2(spikes,ts,handles)
%Calculates the spike features

scales = handles.WC.scales;
feature = handles.WC.features;
exclusioncrit = handles.WC.exclusioncrit; % thr; number
exclusionthr = handles.WC.exclusionthr;
maxinputs = handles.WC.maxinputs;
int_factor = handles.WC.int_factor;
w_pre = handles.WC.w_pre;
w_post = handles.WC.w_post;
% sr = handles.WC.sr;

nspk=size(spikes,1);
len = size(spikes,2)/int_factor;


shift = 2;
offstart = round(w_pre/2) + shift;
offend  = round(w_post/2) + shift + w_pre;

ind1 = [fliplr(w_pre : -1: offstart+1),w_pre+1 : offend] * int_factor;
% ind1 = 8:23
% downfac = 2;
% ind2 = [fliplr(w_pre : -downfac : offstart+1),w_pre+1 : downfac : offend] * int_factor;

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
    norm_factor_WL=diff(prctile(cc(:),[95,5]));
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
    norm_factor_PCA=diff(prctile(S(:),[95,5]));
    S=S*norm_factor_WL/norm_factor_PCA;
    ccall = cat(2,ccall,S);
    fnall = cat(2,fnall,fn);
    clear fn S
end


%% Raw waveforms

if strfind(feature,'raw')
    %     inddummy = ind1;
    %     inddummy(ind1 == w_pre * int_factor) = [];
    %     spikesdown = spikes(:,inddummy);
    
    spikesdown = spikes(:,ind1);
    
    %features
    fn = cell(1,length(ind1));
    for i=1:length(ind1)
        fn{i}=sprintf('Raw,%d',i);
    end
    
    %     fn(ind1 == w_pre * int_factor) = [];
    
    norm_factor_raw=diff(prctile(spikesdown(:),[95,5]));
    spikesdown=spikesdown*norm_factor_WL/norm_factor_raw;
    ccall = cat(2,ccall,spikesdown);
    fnall = cat(2,fnall,fn);
    clear fn inddummy spikesdown
end

%% time !!!!!!!!!!!!

if strfind(feature,'time')
    norm_factor_t=diff(prctile(ts(:),[95,5]));
    draw=ts*norm_factor_WL/norm_factor_t;
    
    ccall = cat(2,ccall,draw);
    fnall = cat(2,fnall,{'Time'});
end

% integral norm
for k=1:size(ccall,2)
    ccall(:,k)=(ccall(:,k)-mean(ccall(:,k)))/sum(abs(ccall(:,k)-mean(ccall(:,k))))*size(ccall,1); % not sure about precision problems here
end



%% feature selection

warning('off','stats:lillietest:OutOfRangeP');
warning('off','stats:lillietest:OutOfRangePHigh');
warning('off','stats:lillietest:OutOfRangePLow');

% stdcc = std(ccall) * 3;
% meancc = mean(ccall);
% thr_dist_min = meancc - stdcc;
% thr_dist_max = meancc + stdcc;
% 
% sd=nan(1,size(ccall,2));
% for i=1:size(ccall,2)                            % lilliefors test for coefficient selection
%     aux = ccall(ccall(:,i) > thr_dist_min(i) & ccall(:,i) < thr_dist_max(i),i); % removing outliers (3x std)
%     if length(aux) > 10;
%         [~,~,lstat]=lillietest(aux,0.05,'ev');
%         sd(i)=lstat;
%     else
%         sd(i)=0;
%     end
% end
% if strfind(feature,'time')
% sd(end)=max(sd)+1; %% to alway keep time??
% end
% 
% [sortsd,ind]=sort(-sd); % orders indexes by lstat (KSstatistic) of lilliefors test, the higher ones first



ccall_reduced=ccall;
n_features=size(ccall_reduced,2);
feature_index=1:n_features;
removed_features=[];
sortsd=[];
while n_features>0
    lstat=nan(1,n_features);
    for f=1:n_features
        
        ccall_tmp=ccall_reduced;
        if n_features>1
        ccall_tmp(:,f)=[];
        end
        mean_feat=mean(ccall_tmp,1);
       [ mode_feat,M]=mode(round(ccall_tmp*100),1);
        distances_multidim=bsxfun(@minus, ccall_tmp, mode_feat);
        %distances_euclidean=sqrt(sum(distances_multidim.^2,2));
        distances_euclidean=mean(distances_multidim,2);
        
%     aux = ccall(ccall(:,i) > thr_dist_min(i) & ccall(:,i) < thr_dist_max(i),i); % removing outliers (3x std)
        aux=distances_euclidean;                %% remove outliers still ?? how?
        [~,~,lstat(f)]=lillietest(aux,0.05,'ev');
    end
    [sd, to_remove]=max(lstat);
    
    %% dont remove time!
    if strfind(feature,'time') && to_remove == size(ccall_reduced,2) && n_features>1
       [sd, to_remove]=max(lstat(1:end-1)); 
    end
    ccall_reduced(:,to_remove)=[]; %try doing it in one dimension only for similar sds
    
    F=feature_index(to_remove);
    removed_features=[removed_features F];
    sortsd=[sortsd sd];
    feature_index(to_remove)=[];
    n_features=numel(feature_index);
end
ind=fliplr(removed_features);


if strfind(feature,'time')
    t_idx=size(ccall,2);
    ind(ind==t_idx)=[];
    ind=[t_idx,ind];
end

%[sortsd,ind]=sort(-sd); % orders indexes by lstat (KSstatistic) of lilliefors test, the higher ones first

















% 
% 
% 
% 
% 
% % excluding features if KStatistics below certain threshold (mean +
% % 1.96*std) of a bootstrapped random dummie
% if strcmp(exclusioncrit,'thr')
%     crit = nan(100,1);
%     for i = 1 : 100
%         dummy = randn(size(spikes,1),1);
%         dummy = dummy(dummy > -3 & dummy < 3);
%         [~,~,lstat]=lillietest(dummy);
%         crit(i) = lstat;
%     end
%     crit = mean(crit) + 1.96*std(crit);
%     
%     ind = ind(sortsd < -crit);
%     %     ind = ind2(sortsd < -0.05);
%     
%     %     if length(ind) > maxinputs
%     %         ind = ind(1 : maxinputs);
%     %     end
% elseif strcmp(exclusioncrit,'sta')
%    
%     %% threshold
%     [sortsd,ind]=sort(-sd); % orders indexes by lstat (KSstatistic) of lilliefors test, the higher ones first
%     
%     crit = nan(100,1);
%     for i = 1 : 100
%         dummy = randn(size(spikes,1),1);
%         dummy = dummy(dummy > -3 & dummy < 3);
%         [~,~,lstat]=lillietest(dummy);
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
% 
% %% removing correlated features! we should probably be more strict here -> lower exclusionthr
% if ~isempty(ccall)
%     feat_corr = corr(ccall).^2;                    % feature correlation matrix
%     feat_corr_true = feat_corr >= exclusionthr;           % logical matrix of correlated features
%     feat_corr_true(diag(true(1,size(feat_corr_true,1)))) = 0; % removing autocorrelation
%     ind = 1 : size(feat_corr_true,1);
%     indrej = [];
%     
%     %% removing features one by one, starting with the one that correlates with most
%     while sum(feat_corr_true(:))
%         distdummy = sum(feat_corr_true);
%         %thrdummy = exclusionthr;
%         delind = find(distdummy == max(distdummy));
%         if numel(delind) > 1 % several features with same max correlation
%             simmat3 = feat_corr(delind,:);
%             simmat3(simmat3<exclusionthr)=2;
%             simmat3=round(simmat3*100)/100; % 0.01 stepsize maintained
%             min_above_thr=min(simmat3,[],2);
%             delind(find(min_above_thr==max(min_above_thr),1, 'last'));
%         end
%         %% before: make 0.01 iterations, check each time how many surpass the threshold
%         %     while length(delind) > 1 % several features with same max correlation
%         %         thrdummy = thrdummy+0.01;
%         %         simmat3 = simmat >= thrdummy;
%         %         simmat3(diag(true(1,size(simmat3,1)))) = 0;
%         %         distdummy = sum(simmat3);
%         %         delind = delind(distdummy(delind) == max(distdummy(delind)));
%         %         if ~sum(distdummy(delind))
%         % %             delind = randsample(delind,1);
%         %             delind = delind(end);
%         %         end
%         %     end
%         feat_corr(:,delind) = [];
%         feat_corr(delind,:) = [];
%         feat_corr_true(:,delind) = [];
%         feat_corr_true(delind,:) = [];
%         indrej = cat(2,indrej,ind(delind));
%         ind(delind)= [];
%     end
%     
%     ccall = ccall(:,ind);
%     fnall = fnall(ind);
% end

if size(ccall,2) < maxinputs
    maxinputs = size(ccall,2);
end

inspk = ccall(:,1:maxinputs);
feature_names=fnall(1:maxinputs);
inputs = maxinputs;


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

