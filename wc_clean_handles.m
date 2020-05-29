function handles=wc_clean_handles(handles)
%something is still not cleaned properly, check temperature plot
MAX_CLUS=handles.MAX_CLUS;
handles.plotted(1:MAX_CLUS)=0;
if isfield(handles,'classind'), handles=rmfield(handles,'classind'); end
handles.fixed(1:MAX_CLUS)=0;
set(handles.hfix,'Value',0);
if isfield(handles,'fixed_classind'), handles=rmfield(handles,'fixed_classind'); end
handles.forced(1:MAX_CLUS)=0;
for i=1:MAX_CLUS, handles.classind_unforced{i}=[]; end
handles.rejected=1;
handles.spikes=[];
handles.classtemp=[];
handles.tree=[];
handles.features=[];
handles.feature_names=[];
handles.ts=[];
handles.ts_time=[];

cla(handles.axesTS);
cla(handles.axesAllClusters);
cla(handles.axesClust0);
%cla(handles.axesISI0);
set(handles.hclassify,'Value',0);
set(handles.hclassify,'String','Classify');
%-------
set(handles.htempmatch,'Value',0);
set(handles.htempmatch,'String','Near');
%-------
set(handles.htempmatch2,'Value',0);
set(handles.htempmatch2,'String','Near T');
%-------
set(handles.htempmatch3,'Value',0);
set(handles.htempmatch3,'String','NearT1');


%-------
cla(handles.axesTemp);
if isfield(handles,'hhor'), handles=rmfield(handles,'hhor'); end
if isfield(handles,'hver'), handles=rmfield(handles,'hver'); end
set(handles.htoplot,'Value',1);%temporaraly for PT
for i=1:numel(handles.spikeaxes), 
    cla(handles.spikeaxes(i));
    %cla(handles.isiaxes(i));
    if i<14, set(handles.hclustergroup{i},'Visible','off'); end % cluster 0?
end
set(handles.hsuppl,'Visible','Off');
set(handles.hsuppl2,'Visible','Off');
set(handles.hdetailsfig,'Visible','Off');
handles.plotted(1:MAX_CLUS)=0;
handles.classind={};
handles.fixed(1:MAX_CLUS)=0;
for i=1:MAX_CLUS, handles.fixed_classind{i}=[]; end
handles.forced(1:MAX_CLUS)=0;
for i=1:MAX_CLUS, handles.classind_unforced{i}=[]; end
handles.rejected=1;
handles.spikes=[];
%--------
handles.spikesadd=[];
handles.spikescomp=[];
%--------
handles.classtemp=[];
handles.tree=[];
handles.features=[];
handles.feature_names=[];
handles.ts=[];
handles.ts_time=[];

handles.index=[];
handles.nspk=NaN;
handles.ncl=NaN;
handles.nfeatures=NaN;
handles.min_clus=NaN;
handles.temp=NaN;
handles.sp_time=[];
handles.mean_ss=[];
handles.std_ss=[];
handles.hcol=[];
% 
% if isfield(handles,'hfeatures'),
%     figure(handles.hfeatures);
%     clf reset; 
%     set(handles.hfeatures,'Visible','Off');
% end
