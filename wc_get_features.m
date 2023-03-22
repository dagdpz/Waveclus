function [features,f_names,sd] = wc_get_features(spikes,ts,handles)
%Calculates the spike features
scales      = handles.WC.scales;
feature     = handles.WC.features;
int_factor  = handles.WC.int_factor;
w_pre       = handles.WC.w_pre;
w_post      = handles.WC.w_post;
nspk        = size(spikes,1);
len         = size(spikes,2)/int_factor;

ind1 = 1:(w_post+w_pre);

% CALCULATES FEATURES
features = [];
f_names = {};
norm_factor_WL=1;

%% Wavelet decomposition
if strfind(feature,'wav')
    cc=zeros(nspk,length(ind1));
    for i=1:nspk
        [c,l]=wavedec(spikes(i,ind1),scales,handles.WC.wavelet);
        cc(i,:)=c;
    end
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
    features = cat(2,features,cc);
    f_names = cat(2,f_names,fn);
    norm_factor_WL=std(cc(:));
    clear fn cc
end

%% PCA
if strfind(feature,'pca')
    numdim = 10;
    matlabversion=datevec(version('-date'));
    if matlabversion(1)>=2014
        [~,S] = pca(spikes(:,ind1),'Economy',false); 
    else
        [~,S] = princomp(spikes(:,ind1)); 
    end
    fn = cell(1,numdim);
    for i=1:numdim
        fn{i}=sprintf('PCA,%d',i);
    end
    S=S(:,1:numdim);
    norm_factor_PCA=std(S(:));
    S=S*norm_factor_WL/norm_factor_PCA;
    features = cat(2,features,S);
    f_names = cat(2,f_names,fn);
    clear fn S
end

%% Raw waveforms
if strfind(feature,'raw')
    spikesdown = spikes(:,ind1);
    fn = cell(1,length(ind1));
    for i=1:length(ind1)
        fn{i}=sprintf('Raw,%d',i);
    end
    norm_factor_raw=std(spikesdown(:));
    spikesdown=spikesdown*norm_factor_WL/norm_factor_raw;
    features = cat(2,features,spikesdown);
    f_names = cat(2,f_names,fn);
    clear fn inddummy spikesdown
end

%% time 
if strfind(feature,'time')
    norm_factor_t=1800000; %% normalized to 30 minutes
    draw=ts*norm_factor_WL/norm_factor_t;
    features = cat(2,features,draw);
    f_names = cat(2,f_names,{'Time'});
end

%% feature selection
warning('off','stats:lillietest:OutOfRangeP');
warning('off','stats:lillietest:OutOfRangePHigh');
warning('off','stats:lillietest:OutOfRangePLow');

stdcc = nanstd(features, [], 1) * 3; % calculate std exclusively by row
meancc = nanmean(features, 1); % calculate mean exclusively by row
thr_dist_min = meancc - stdcc;
thr_dist_max = meancc + stdcc;

sd=nan(1,size(features,2));
for i=1:size(features,2)                            % lilliefors test for coefficient selection
    aux = features(features(:,i) > thr_dist_min(i) & features(:,i) < thr_dist_max(i),i); % removing outliers (3x std)
    if length(aux) > 10;
        [~,~,lstat]=lillietest(aux,0.05);
        sd(i)=lstat;
    else
        sd(i)=0;
    end
end

warning('on','stats:lillietest:OutOfRangeP');
warning('on','stats:lillietest:OutOfRangePHigh');
warning('on','stats:lillietest:OutOfRangePLow');
end