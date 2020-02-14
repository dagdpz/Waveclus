function [inspk,feature_names] = wave_features4(spikes,handles)
%Calculates the spike features

scales = handles.par.scales;
feature = handles.par.features;
similarity = handles.par.similarity;
% inputs = handles.par.inputs;
maxinputs = handles.par.maxinputs;
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
    clear fn
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
    clear fn
end


%% Raw waveforms

if strfind(feature,'raw')
    inddummy = ind1;
    inddummy(ind1 == w_pre * int_factor) = [];
    spikesdown = spikes(:,inddummy);
    
    %features
    fn = cell(1,length(ind1));
    for i=1:length(ind1)  %/downfac,
        fn{i}=sprintf('Raw,%d',i);
    end
    
    fn(ind1 == w_pre * int_factor) = [];
    
    ccall = cat(2,ccall,spikesdown);
    fnall = cat(2,fnall,fn);
    clear fn inddummy
end


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

warning('on','stats:lillietest:OutOfRangeP');
warning('on','stats:lillietest:OutOfRangePHigh');
warning('on','stats:lillietest:OutOfRangePLow');

ccall(:,sd < 0.02) = [];
fnall(sd < 0.02) = [];
sd(sd < 0.02) = [];


[~,ind]=sort(-sd);
ccall = ccall(:,ind);
fnall = fnall(ind);

simmat = corr(ccall).^2 >= similarity;
simmat(diag(true(1,size(ccall,2)))) = 0;
% figure;imagesc(simmat)


while sum(simmat(:))
    firstind = find(sum(simmat) > 0,1,'first');
    delind = simmat(:,firstind);
    simmat(:,delind) = [];
    simmat(delind,:) = [];
    ccall(:,delind) =[];
    fnall(delind) =[];
end

if size(ccall,2) > maxinputs
    feature_names={fnall{1:maxinputs}};
    inspk=ccall(:,1:maxinputs);
else
    feature_names={fnall};
    inspk=ccall;
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

