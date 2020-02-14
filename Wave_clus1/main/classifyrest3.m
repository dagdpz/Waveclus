function handles=classifyrest3(handles);
%add possibility to choose between spike shapes and features 

handles.par.template_sdnum = 5;             % max radius of cluster in std devs.
shift = 2;
shift = shift*handles.par.int_factor;

classes=zeros(1,handles.nspk);
for i=1:handles.ncl
    classes(handles.classind{i})=i;
end
if ~sum(classes==0), return; end %nothing to do, maybe check for this earlier
group=classes([classes~=0]);
class0=find(classes==0);
switch handles.classify_space,
    case 'spikeshapes', 
        ints=handles.par.w_pre*handles.par.int_factor-handles.par.w_pre*handles.par.int_factor/2+1:handles.par.w_pre*handles.par.int_factor+handles.par.w_post*handles.par.int_factor/2;
        sample=handles.spikes(class0,ints);
        training=handles.spikes(classes~=0,ints);
    case 'features',
        sample=handles.features(class0,:);
        training=handles.features(classes~=0,:);
    case 'spikeshapesfeatures',
        ints1 = sort(handles.par.w_pre*handles.par.int_factor: -handles.par.int_factor: handles.par.w_pre*handles.par.int_factor-handles.par.w_pre*handles.par.int_factor/2+1+shift);
        ints2 = handles.par.w_pre*handles.par.int_factor + handles.par.int_factor : handles.par.int_factor: handles.par.w_pre*handles.par.int_factor+handles.par.w_post*handles.par.int_factor/2+shift;
        ints = cat(2,ints1,ints2);
        clear ints1 ints2
%         sample=cat(2,handles.spikes(class0,ints),handles.spikescomp(class0,:),handles.spikesadd(class0,:));
%         training=cat(2,handles.spikes(classes~=0,ints),handles.spikescomp(classes~=0,:),handles.spikesadd(classes~=0,:));
%         sample=cat(2,handles.spikes(class0,ints),handles.spikesadd(class0,:));
%         training=cat(2,handles.spikes(classes~=0,ints),handles.spikesadd(classes~=0,:));
        sample=cat(2,handles.spikes(class0,ints));
        training=cat(2,handles.spikes(classes~=0,ints));
        
    case 'selectedfeatures',
        
end
switch handles.classify_method,
    case 'force',
        %rodrigo's methods
    otherwise
        
        [centers, sd, pd] = build_templates(group,training); % we are going to ignore pd
        sdnum = handles.par.template_sdnum;
        
        index = false(length(class0),1);
        for i = 1 : length(class0)
            distances = sqrt(sum((ones(size(centers,1),1)*sample(i,:) - centers).^2,2)');
            conforming = find(distances < sdnum*sd);
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
            [c,err] = classify(sample, training, group, handles.classify_method);
        else
            c=ones(1,size(sample,1));
        end
        c(index) = 0;
end
handles.classind_unforced=handles.classind;
for i=1:handles.ncl,
    handles.classind{i}=sort([handles.classind{i} class0(c==i)]);
    handles.forced(i)=1;
end
handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);
