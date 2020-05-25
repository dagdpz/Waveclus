function outclass = wc_classify2(sample, training, group)
carry_on=true;
outidxes=1:size(sample,1);
unique_groups=unique(group);
outclass=zeros(size(outidxes));
while carry_on
    carry_on=false;
    group_prevalence=hist(group,unique_groups);
    samples_to_classify=[];
    groups_to_classify=[];
    for g=unique(group)
        stepsize_per_group(g,:)=std(training(group==g,:),1,1);
    end
    stepsize=stepsize_per_group(group,:);
    for s=1:size(sample,1)
        distance=abs(bsxfun(@minus,training,sample(s,:)));
        distance=distance-stepsize;
        within_distance=distance<0;
        idx_close=all(within_distance,2);
        if any(idx_close)
            close_groups=group(idx_close);
            hh=hist(close_groups,unique_groups)./group_prevalence;
            [~,g]=max(hh);
            g=unique_groups(g);
            samples_to_classify=[samples_to_classify;s];
            groups_to_classify=[groups_to_classify;g];
            carry_on=true;
        end
    end
    outclass(outidxes(samples_to_classify))=groups_to_classify;
    training=[training;sample(samples_to_classify,:)];
    group=[group groups_to_classify'];
    outidxes(samples_to_classify)=[];
    sample(samples_to_classify,:)=[];
end
end

