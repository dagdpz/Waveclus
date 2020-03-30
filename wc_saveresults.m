function handles=wc_saveresults(handles)
par=handles.par;
par.temp=handles.temp;
par.min_clus=handles.min_clus;

cluster_class=zeros(handles.nspk,2);
cluster_class(:,2)= handles.index;
for i=1:handles.ncl, cluster_class(handles.classind{i},1)=i; end
%save data
timesfile=sprintf('%s.mat',handles.filename);
save(timesfile,'par','cluster_class','-append');
%save pictures
name=sprintf('ch%d-%s.jpg',handles.channel,handles.bname);
name1=sprintf('ch%d-%s_b.jpg',handles.channel,handles.bname);
name2=sprintf('ch%d-%s_c.jpg',handles.channel,handles.bname);

if strcmp(get(handles.hsuppl,'Visible'),'on'), 
    name=sprintf('ch%d-%s_a.jpg',handles.channel,handles.bname);
    figure(handles.hsuppl);
    print(gcf,'-djpeg',name1);
end
if strcmp(get(handles.hsuppl2,'Visible'),'on'), 
    figure(handles.hsuppl2);
    print(gcf,'-djpeg',name2);
end
print(handles.mainfig,'-djpeg',name);

if isfield(handles,'hfeatures'),
    if strcmp(get(handles.hfeatures,'Visible'),'on'), 
        figure(handles.hfeatures);
        name=sprintf('ch%d-%s-features.jpg',handles.channel,handles.bname);
        print(gcf,'-djpeg',name);
    end
end
