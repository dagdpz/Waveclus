function handles=tempmatch(handles);
%add possibility to choose between spike shapes and features 
handles.par.template_type = 'center';       % nn, center, ml, mahal
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

% f_in  = handles.spikes(sort(cell2mat(handles.classind(1,1:end-1))),:);
% f_out = handles.spikes(cell2mat(handles.classind(1,end)),:);

% f_in  = handles.features(sort(cell2mat(handles.classind(1,1:end-1))),:);
% f_out = handles.features(cell2mat(handles.classind(1,end)),:);
% 
% f_in  = cat(2,handles.spikes(sort(cell2mat(handles.classind(1,1:end-1))),preoff+1:end-postoff),handles.features(sort(cell2mat(handles.classind(1,1:end-1))),1:postoffB));
% f_out = cat(2,handles.spikes(cell2mat(handles.classind(1,end)),preoff+1:end-postoff),handles.features(cell2mat(handles.classind(1,end)),1:postoffB));


f_in  = cat(2,handles.spikes(sort(cell2mat(handles.classind(1,1:end-1))),handles.par.w_pre*handles.par.int_factor-handles.par.w_pre*handles.par.int_factor/2+1+shift:handles.par.w_pre*handles.par.int_factor+handles.par.w_post*handles.par.int_factor/2+shift),handles.spikescomp(sort(cell2mat(handles.classind(1,1:end-1))),:),handles.spikesadd(sort(cell2mat(handles.classind(1,1:end-1))),:));
f_out = cat(2,handles.spikes(cell2mat(handles.classind(1,end)),handles.par.w_pre*handles.par.int_factor-handles.par.w_pre*handles.par.int_factor/2+1+shift:handles.par.w_pre*handles.par.int_factor+handles.par.w_post*handles.par.int_factor/2+shift),handles.spikescomp(cell2mat(handles.classind(1,end)),:),handles.spikesadd(cell2mat(handles.classind(1,end)),:));


c = force_membership_wc(f_in, group, f_out,handles);

handles.classind_unforced=handles.classind;
for i=1:handles.ncl,
    handles.classind{i}=sort([handles.classind{i} class0(c==i)]);
    handles.forced(i)=1;
end
handles.classind{end}=setdiff(1:handles.nspk,[handles.classind{1:end-1}]);