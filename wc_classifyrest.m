function handles=wc_classifyrest(handles,method)
%add possibility to choose between spike shapes and features
%rework entirely......
shift = 2;
shift = shift*handles.WC.int_factor;

classes=zeros(1,handles.nspk);
for i=1:handles.ncl
    classes(handles.classind{i})=i;
end
if ~sum(classes==0), return; end %nothing to do, maybe check for this earlier
group=classes([classes~=0]);
class0=find(classes==0);
switch handles.WC.classify_space,
    case 'spikeshapes',
        ints=handles.WC.w_pre*handles.WC.int_factor-handles.WC.w_pre*handles.WC.int_factor/2+1:handles.WC.w_pre*handles.WC.int_factor+handles.WC.w_post*handles.WC.int_factor/2;
        sample=handles.spikes(class0,ints);
        training=handles.spikes(classes~=0,ints);
    case 'features',
        sample=handles.features(class0,:);
        training=handles.features(classes~=0,:);
    case 'spikeshapesfeatures',
        ints1 = sort(handles.WC.w_pre*handles.WC.int_factor: -handles.WC.int_factor: handles.WC.w_pre*handles.WC.int_factor-handles.WC.w_pre*handles.WC.int_factor/2+1+shift);
        ints2 =      handles.WC.w_pre*handles.WC.int_factor + handles.WC.int_factor : handles.WC.int_factor: handles.WC.w_pre*handles.WC.int_factor+handles.WC.w_post*handles.WC.int_factor/2+shift;
        ints = cat(2,ints1,ints2);
        clear ints1 ints2
        sample=cat(2,handles.spikes(class0,ints));
        training=cat(2,handles.spikes(classes~=0,ints));
    case 'selectedfeatures',
        
end
switch handles.WC.classify_method,
    case 'force',
        %rodrigo's methods
    otherwise
        [centers, sd] = wc_build_templates(group,training);
        sdnum = handles.WC.template_sdnum;
        index = false(length(class0),1);
        for i = 1 : length(class0) % each spike
            distances = sqrt(sum((ones(size(centers,1),1)*sample(i,:) - centers).^2,2)'); % check distance in all features
            conforming = distances < sdnum*sd;
            if sum(conforming)==0 % spike is off in all features
                index(i) = 1;     % will be put in 0 cluster in the end
            end
            clear distances conforming
        end
        if handles.ncl>1,
            while size(training,2)>size(training,1),%if number of features is bigger then number of spikes % when???
                sample=sample(:,1:2:end);
                training=training(:,1:2:end);
            end
            if method==1
                c = wc_classify2(sample, training,group);
            elseif method==2
                c = wc_classify(sample, training,group, handles.WC.classify_method);
            end
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
