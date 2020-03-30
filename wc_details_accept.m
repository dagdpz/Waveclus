function handles=wc_details_accept(handles)
%remove bad spikes if neccesary, put them into cluster 0
%fix cluster?

handles.classind{end}=unique([setdiff(handles.classind{handles.details_currcluster},handles.details_inds_accepted) handles.classind{end}]);
handles.classind{handles.details_currcluster}=sort(handles.details_inds_accepted);
