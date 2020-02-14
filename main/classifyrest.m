function handles=classifyrest(handles);
%add possibility to choose between spike shapes and features 

classes=zeros(1,handles.nspk);
for i=1:handles.ncl
    classes(handles.classind{i})=i;
end
if ~sum(classes==0), return; end %nothing to do, maybe check for this earlier
group=classes([classes~=0]);
class0=find(classes==0);
switch handles.classify_space,
    case 'spikeshapes', 
        ints=handles.par.w_pre-handles.par.w_pre/2:handles.par.w_pre+handles.par.w_post/2-1;
        sample=handles.spikes(class0,ints);
        training=handles.spikes(classes~=0,ints);
    case 'features',
        sample=handles.features(class0,:);
        training=handles.features(classes~=0,:);
    case 'selectedfeatures',
        
end
switch handles.classify_method,
    case 'force',
        %rodrigo's methods
    otherwise
        if handles.ncl>1,
            while size(training,2)>size(training,1),%if number of features is bigger then number of spikes
                sample=sample(:,1:2:end);
                training=training(:,1:2:end);
            end
            [c,err] = classify(sample, training, group, handles.classify_method);
        else
            c=ones(1,size(sample,1));
        end
end
handles.classind_unforced=handles.classind;
for i=1:handles.ncl,
    handles.classind{i}=sort([handles.classind{i} class0(c==i)]);
    handles.forced(i)=1;
end
handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);
