function handles=get_clustermean(handles);
%estimates a cluster mean in predifined proximity space
%also estimates the distances from each spike in a cluster to the mean

i=handles.details_currcluster;
%do some preprocessing for current cluster
%find center of cluster in feature or spike shape space
if handles.forced(i), inds=handles.classind_unforced{i};
else inds=handles.classind{i}; end

handles.details_mean_spikesshape=mean(handles.spikes(inds,:),1);


switch handles.classify_space,
    case 'features',%all features weighted equally
        %if cluster was forced (aka template matching was used) then use an
        %unforced version to estimate the center
        if handles.forced(i), inds=handles.classind_unforced{i};
        else inds=handles.classind{i}; end

        %estimate center
        mt=mean(handles.features(inds,:),1);
        %estimate distance to the center
        inds=handles.classind{i};
        t=handles.features(inds,:);
        st=[length(inds) 1];
        t=(t-repmat(mt,st))./repmat(std(t,[],1),st);
        handles.details_dist2center=t;
    case 'selectedfeatures',%NOT DONE%weight features according to their importance for current cluster
        %find the most separating features
        if handles.forced(i), inds=handles.classind_unforced{i};
        else inds=handles.classind{i}; end
        inds0=handles.classind{end};
        [nspk,nf]=size(handles.features);
        notinds=setdiff(1:nspk,[inds inds0]);
        if length(notinds)>1,
            for j=1:nf,
                zval(j)=NaN;
                [prs(j),h,stats]=ranksum(handles.features(inds,j),handles.features(notinds,j));
                if isfield(stats,'zval'), zval(j)=-abs(stats.zval); end
            end
            [t,indrs]=sort(zval);
        end
        %normalize featues
        t=handles.features(handles.classind{i},:);
        st=[size(t,1) 1];
        mt=mean(t,1);
        t=(t-repmat(mt,st))./repmat(std(t,[],1),st);
        handles.details_dist2center=t;
    case 'spikeshapes',
        %take whole spike shape
        %         mt=mean(handles.spikes(handles.classind{i},:),1);
        %         handles.details_dist2center=handles.spikes(handles.classind{i},:)-repmat(mt,[length(handles.classind{i}) 1]);
        
        %take only half around min(max)
        ints=handles.par.w_pre-handles.par.w_pre/2:handles.par.w_pre+handles.par.w_post/2-1;
        %         mt=mean(handles.spikes(handles.classind{i},ints),1);
        mt=mean(handles.spikes(handles.classind_unforced{i},ints),1);
        handles.details_dist2center=handles.spikes(handles.classind{i},ints)-repmat(mt,[length(handles.classind{i}) 1]);

end
%square distance to the center of cluster
handles.details_dist2center=sum(handles.details_dist2center.^2,2);
[t,ind]=sort(handles.details_dist2center);
%indices of spikes sorted according to the distance to the center of cluster 
handles.details_inds_ofdistsorted=handles.classind{i}(ind);
