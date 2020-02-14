function combine_clusters(filename, combs)

if ~iscell(combs),
    t=combs;
    clear combs
    combs{1}=t;
end
q=load(filename);
cc=q.cluster_class(:,1);
for ic=1:length(combs),
    c=combs{ic};
    for j=2:length(c),
        cc(q.cluster_class(:,1)==c(j))=c(1);
    end
end
uc=unique(cc);
if uc(1)==0, uc=uc(2:end); end
cc1=cc;
for i=1:length(uc),
    if i==uc(i), continue; end
    cc1(cc==uc(i),1)=i;
end
cluster_class=q.cluster_class;
cluster_class(:,1)=cc1;
save(filename,'cluster_class','-append')
    
