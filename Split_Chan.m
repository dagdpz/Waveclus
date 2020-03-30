function Split_Chan(handles)

dtype = handles.dtyperead;
numchan = handles.numchan;
rawname = handles.rawname;


files = dir(rawname);
filesB = files.name;

fid=fopen(filesB,'r','ieee-le');
dataall = fread(fid,inf,[dtype '=>' dtype]);
fclose(fid); 
datalen = length(dataall)/numchan;
dataall = reshape(dataall,[numchan datalen]);


for i = 1 : numchan
    data = dataall(i,:)';
    save(['dataraw_ch' sprintf('%03d',i)],'data')
    clear data
end


