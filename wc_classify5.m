function outclass = wc_classify5(sample, training, group)
carry_on=true;
outidxes=1:size(sample,1);
unique_groups=unique(group);
outclass=zeros(size(outidxes));
features_idx=1:size(training,2);
while carry_on
    carry_on=false;
    group_prevalence=hist(group,unique_groups);
    samples_to_classify=[];
    groups_to_classify=[];
    u_groups=unique(group);
    
    for g=u_groups
        mean_per_group(g,:)=mean(training(group==g,:),1);
    end
    for g=u_groups
        sep_per_group(g,:)=min(abs(bsxfun(@minus,mean_per_group(u_groups~=g,:),mean_per_group(g,:))));
    end
    for g=u_groups
        stepsize_per_group(g,:)=std(training(group==g,:),0,1);
        %% new
        [~,min_idx]=sort(-sep_per_group(g,2:end));
        max_idx=features_idx(~ismember(features_idx,[1 min_idx(1:4)+1])); %% take only 4 dimensions?
        stepsize_per_group(g,max_idx)=Inf;
    end
    stepsize=stepsize_per_group(group,:);
    UB=training+stepsize;
    LB=training-stepsize;
    within=false(size(training,1),size(sample,1));
    hh=zeros(size(unique_groups,2),size(sample,1));
    groups_within=zeros(size(training,1),size(sample,1));
    for s=1:size(sample,1)
       within(:,s)= all(bsxfun(@ge,UB,sample(s,:)) & bsxfun(@le,LB,sample(s,:)),2);
       groups_within( within(:,s),s)=group( within(:,s));
       hh(:,s)=hist(group( within(:,s)),unique_groups)./group_prevalence;
    end
    
    idx=any(hh,1);
    [~,g]=max(hh,[],1);
    if any(idx)
        samples_to_classify=find(idx);
        groups_to_classify=g(idx);
        
        
        outclass(outidxes(samples_to_classify))=groups_to_classify;
        training=[training;sample(samples_to_classify,:)];
        group=[group groups_to_classify];
        outidxes(samples_to_classify)=[];
        sample(samples_to_classify,:)=[];
        
        carry_on=true;
    end
            
%     %% this one
%             stepsize_sample=std(sample,1,1);
%     for s=1:size(sample,1)
%         distance=abs(bsxfun(@minus,training,sample(s,:)));
%         distance=distance-stepsize;
%         within_distance=distance<0;
%         idx_close=all(within_distance,2);
%         if any(idx_close)
%             %% 0 cluster competing
% %             distance_sample=abs(bsxfun(@minus,sample,sample(s,:)));
% %             distance_sample(s,:)=[];
% %             distance_sample=distance_sample-repmat(stepsize_sample,size(distance_sample,1),1);
% %             unclustered_prev=sum(all(distance_sample<0,2))/numel(sample);
% %             
%             close_groups=group(idx_close);
%             hh=hist(close_groups,unique_groups)./group_prevalence;
%             [max_prev,g]=max(hh);
%         %    if max_prev>unclustered_prev
%             samples_to_classify=[samples_to_classify;s];
%             groups_to_classify=[groups_to_classify;g];
%             carry_on=true;
%          %   end
%         end
%     end
end
end

