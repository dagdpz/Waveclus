function handles=classifyrestoverlayed2(handles);
%add possibility to choose between spike shapes and features 

[centers, sd, pd] = build_templates(handles.group,handles.training); % we are going to ignore pd

index = false(size(handles.sample,1),1);
for i = 1 : size(handles.sample,1)
    distances = sqrt(sum((ones(size(centers,1),1)*handles.sample(i,:) - centers).^2,2)');
    conforming = find(distances < handles.sdnum*sd);
    if isempty(conforming)
        index(i) = 1;
    end
    clear distances conforming
end

if handles.ncl>1,
    while size(handles.training,2)>size(handles.training,1),%if number of features is bigger then number of spikes
        handles.sample=handles.sample(:,1:2:end);
        handles.training=handles.training(:,1:2:end);
    end
    [handles.c,err] = classify(handles.sample, handles.training, handles.group, handles.classify_method);
else
    handles.c=ones(1,size(handles.sample,1));
end
handles.c(index) = 0;


