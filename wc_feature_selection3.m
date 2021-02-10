function [inspk,feature_names,inputs] = wc_feature_selection3(ccall,fnall,feature_sds,features_per_subset,handles)
%% get relevant parameters
feature = handles.WC.features;
exclusioncrit = handles.WC.exclusioncrit; % thr; number
exclusionthr = handles.WC.exclusionthr;
maxinputs = handles.WC.maxinputs;

%% order features so that top features within subset become top features over all 
feats_to_order=fnall;
idx=0;
while numel(feats_to_order)>0
    for k=1:size(feature_sds,1)
        idx=idx+1;
        if numel(feats_to_order)==0
            continue;
        end
        [~,f]=max(feature_sds(k,:));        
        ordered_features{idx}=feats_to_order{f};
        feature_sds(:,f)=[];
        feats_to_order(f)=[];
    end
end


%% compute average features of each spike subset (SUs, MUs, both - potentially with neg and pos threshold)
for k=1:numel(features_per_subset)    
    av_features(k,:)=mean(features_per_subset{k},1);
end    
if numel(features_per_subset)<=2 %% only one threshold
    compare_avarages=[1 2];
else %% two thresholds
    compare_avarages=[1 2; 3 4; 1 3];
end

top_features={'Time'};%% to always keep time
% if strfind(feature,'pca') %% to always keep pca 1,2,3 
%     top_features=[top_features {'PCA,1','PCA,2','PCA,3'}];
% end
for comp=1:size(compare_avarages,1)
    [~,idx]=max(abs(av_features(compare_avarages(comp,1),:)-av_features(compare_avarages(comp,2),:)));
    top_features{end+1}=fnall{idx};
end

%% put top_features at the beginning!
ordered_features=[top_features ordered_features(~ismember(ordered_features,top_features))];
[~,ind]=ismember(ordered_features,fnall);

ccall = ccall(:,ind);
fnall = fnall(ind);

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



clear ind

%% removing correlated features! we should probably be more strict here -> lower exclusionthr
if ~isempty(ccall)
    feat_corr = corr(ccall).^2;                    % feature correlation matrix
    %% ignore top features!
    %feat_corr(1:numel(top_features),:)=deal(0);   
    
    feat_corr_true = feat_corr >= exclusionthr;           % logical matrix of correlated features
    
    [A,B]=find(feat_corr_true);
    
    f_to_remove=[];
    for cc=unique(A)'
        if ismember(cc, f_to_remove)
            continue;
        end
            
        possible=B(A==cc);
        f_to_remove=[f_to_remove; possible(possible>cc)];
    end
    
    f_to_remove=unique(f_to_remove);
    ccall(:,f_to_remove)=[];
    fnall(f_to_remove)=[];

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
