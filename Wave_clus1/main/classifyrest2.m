function handles=classifyrest2(handles);
%add possibility to choose between spike shapes and features 

handles.par.template_sdnum = 4;             % max radius of cluster in std devs.


classes=zeros(1,handles.nspk);
for i=1:handles.ncl
    classes(handles.classind{i})=i;
end
if ~sum(classes==0), return; end %nothing to do, maybe check for this earlier
group=classes([classes~=0]);
class0=find(classes==0);
switch handles.classify_space,
    case 'spikeshapes', 
        ints=handles.par.w_pre-handles.par.w_pre/2+1:handles.par.w_pre+handles.par.w_post/2;
        sample=handles.spikes(class0,ints);
        training=handles.spikes(classes~=0,ints);
    case 'features',
        sample=handles.features(class0,:);
        training=handles.features(classes~=0,:);
    case 'spikeshapesfeatures', 
        ints=round(handles.par.w_pre-handles.par.w_pre/2:handles.par.w_pre+handles.par.w_post/2-1);
        ints2 = [];
        for i = 1 : handles.par.scales - 1 
            divfac = 2^(i + 1);
            divfac2 = 2^i;
            offfac = size(handles.spikecomp,2)/(2^i);
            intsB = round((handles.par.w_pre/divfac2-handles.par.w_pre/divfac:handles.par.w_pre/divfac2+handles.par.w_post/divfac-1)+offfac);
            ints2 = cat(2,ints2,intsB);
        end
        intsB = round(handles.par.w_pre/divfac2-handles.par.w_pre/divfac:handles.par.w_pre/divfac2+handles.par.w_post/divfac-1);
        ints2 = cat(2,ints2,intsB);
        ints2 = sort(ints2);
%         badfeatures = false(length(handles.par.inputs),1);
%         for i = 1 : handles.par.inputs
%             badfeatures(i) = ~isempty(strfind(handles.feature_names{1,i},'S'));
%         end
        
        sample=cat(2,handles.spikes(class0,ints),handles.spikecomp(class0,ints2));
        training=cat(2,handles.spikes(classes~=0,ints),handles.spikecomp(classes~=0,ints2));
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
