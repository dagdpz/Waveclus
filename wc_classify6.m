function outclass = wc_classify6(sample, training, group)
carry_on=true;
outidxes=1:size(sample,1);
unique_groups=unique(group);
outclass=zeros(size(outidxes));
features_idx=1:size(training,2);
while carry_on
    carry_on=false;
    std_per_feature=std([training],1);
    std_per_feature(1)=std_per_feature(1)/10; %% making time mor eimportant
    group_prevalence=hist(group,unique_groups);
    classified=[];
    histo=zeros([size(sample,1) numel(unique_groups)]);
    for s=1:size(sample,1)
        std_thr=1/2*numel(features_idx);
        N_above_max=2;
        while N_above_max>1
            within= bsxfun(@minus,training,sample(s,:)).^2*(1./std_per_feature'.^2)<std_thr;
            groups_within=group( within);
            hh=hist(groups_within,unique_groups)./group_prevalence;
            N_above_max=sum(hh>1/10)+sum(hh==max(hh) & hh>0)-1;
            std_thr=std_thr*0.5;
            histo(s,:)=hh;
        end
%         if any(hh)
%             carry_on=true;
%             [~,g]=max(hh./group_prevalence);
%             outclass(outidxes(s))=g;
%             training=[training;sample(s,:)];
%             group=[group g];
%             classified=[classified s];
%             
%         end
        
    end
    
    
    idx=any(histo,2);
    [~,g]=max(histo,[],2);
    if any(idx)
        samples_to_classify=find(idx);
        groups_to_classify=g(idx)';
        
        
        outclass(outidxes(samples_to_classify))=groups_to_classify;
        training=[training;sample(samples_to_classify,:)];
        group=[group groups_to_classify];
        outidxes(samples_to_classify)=[];
        sample(samples_to_classify,:)=[];
        
        carry_on=true;
     end
%     
%     sample(classified,:)=[];
%     outidxes(classified)=[];
    
    
end
end

