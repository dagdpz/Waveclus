function outclass = wc_classify4(sample, training, group)
carry_on=true;
outidxes=1:size(sample,1);
unique_groups=unique(group);
outclass=zeros(size(outidxes));
for g=unique(group)
    stepsize_per_group(g,:)=std(training(group==g,:),1,1);
    template_per_group(g,:)=mean(training(group==g,:),1);
end
    samples_to_classify=[];
    groups_to_classify=[];
for s=1:size(sample,1)
    distance=abs(bsxfun(@minus,template_per_group,sample(s,:)));
    distance=sqrt(sum(distance.^2,2));
    [~,g]=min(distance);
    samples_to_classify=[samples_to_classify;s];
    groups_to_classify=[groups_to_classify;g];
end
outclass(outidxes(samples_to_classify))=groups_to_classify;
training=[training;sample(samples_to_classify,:)];
group=[group groups_to_classify'];
outidxes(samples_to_classify)=[];
sample(samples_to_classify,:)=[];
end

