function handles=wc_saveresults(handles)
par=handles.WC;
par.temp=handles.temp;
par.min_clus=handles.min_clus;

cluster_class=zeros(handles.nspk,2);
cluster_class(:,2)= handles.index;
for i=1:handles.ncl, cluster_class(handles.classind{i},1)=i; end
%save data
timesfile=sprintf('%s.mat',handles.filename);
save([handles.pathname filesep timesfile],'par','cluster_class','-append');
%save pictures
indx=strfind(handles.filename,'_ch');
fname=['Ch-' handles.filename(indx+3:end)];
name=sprintf('%s-clusters.jpg',fname);
name1=sprintf('%s-clustersb.jpg',fname);
name2=sprintf('%s-clustersc.jpg',fname);

print(handles.mainfig,'-djpeg',[handles.pathname filesep name]);
if strcmp(get(handles.hsuppl,'Visible'),'on'), 
    figure(handles.hsuppl);
        name=[handles.pathname filesep name1];
    print(gcf,'-djpeg',name);
end
if strcmp(get(handles.hsuppl2,'Visible'),'on'), 
    figure(handles.hsuppl2);
        name=[handles.pathname filesep name2];
    print(gcf,'-djpeg',name);
end

if isfield(handles,'hfeatures'),
    if strcmp(get(handles.hfeatures,'Visible'),'on'), 
        figure(handles.hfeatures);
        name=[handles.pathname filesep fname '-features'];
        print(gcf,'-djpeg',name);
    end
end
if isfield(handles,'htimecourse'),
    if strcmp(get(handles.htimecourse,'Visible'),'on'), 
        figure(handles.htimecourse);
        name=[handles.pathname filesep fname '-timecourse'];
        print(gcf,'-djpeg',name);
    end
end
