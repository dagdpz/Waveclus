function handles=classifyrestoverlayed3(handles,training,sample,group);
%add possibility to choose between spike shapes and features 

[centers, sd, pd] = build_templates(group,training); % we are going to ignore pd

index = false(size(sample,1),1);
for i = 1 : size(sample,1)
    distances = sqrt(sum((ones(size(centers,1),1)*sample(i,:) - centers).^2,2)');
    conforming = find(distances < handles.sdnum*sd);
    if isempty(conforming)
        index(i) = 1;
    end
    clear distances conforming
end

if handles.ncl>1,
    while size(training,2)>size(training,1),%if number of features is bigger then number of spikes
        sample=sample(:,1:2:end);
        training=training(:,1:2:end);
    end
    [handles.c,err] = classify(sample, training, group, handles.classify_method);
else
    handles.c=ones(1,size(sample,1));
end
handles.c(index) = 0;

