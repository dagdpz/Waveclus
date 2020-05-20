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
    for s=1:size(sample,1)
        std_thr=1/2*numel(features_idx);
        N_above_max=2;
        while N_above_max>1
            within= bsxfun(@minus,training,sample(s,:)).^2*(1./std_per_feature'.^2)<std_thr;
            groups_within=group( within);
            hh=hist(groups_within,unique_groups);
            N_above_max=sum(hh./group_prevalence>1/10)+sum(hh==max(hh) & hh>0)-1;
            std_thr=std_thr*0.5;
        end
        if any(hh)
            carry_on=true;
            [~,g]=max(hh./group_prevalence);
            outclass(outidxes(s))=g;
            training=[training;sample(s,:)];
            group=[group g];
            classified=[classified s];
            
        end
        
    end
    sample(classified,:)=[];
    outidxes(classified)=[];
    
    
end
end

