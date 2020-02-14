function SpikedelArt8B(badchan,handles)

numchan = handles.numchan;
numArray = handles.numArray;

chunks = 1000000;
thresholdtodel = 0.36;

warning off

chanindall = 1 : numchan * numArray;
chanindall = setdiff(chanindall,badchan);
datalen = whos('-file',['datafilt_ch' sprintf('%03d',chanindall(1))],'-regexp','data');
datalen = max(datalen.size);

% % ------------ old data --------------
% fid_read = fopen(['datahpch' sprintf('%03d',chanindall(1))]);
% fseek(fid_read, 0, 'eof');
% datalen = ftell(fid_read)/2;  
% fclose(fid_read);
% % ----------------------------------

fid_dummy = nan(1,length(chanindall));
for i = 1 : length(chanindall)
    fid_dummy(i) = fopen(['dummy' sprintf('%03d',chanindall(i))],'w+');
end

dummyoffset = 0;

for j = 1 : numArray
    tic
    
    chanind = chanindall(chanindall > (j-1) * numchan & chanindall <= j * numchan);
    
    dataall=zeros(datalen,length(chanind),'single');
    for i = 1 : length(chanind)
        load(['datafilt_ch' sprintf('%03d',chanind(i))],'data');
        
%         % ------------ old data --------------
%         fid_read = fopen(['datahpch' sprintf('%03d',chanind(i))]);
%         data = fread(fid_read,'int16=>int16');
%         fclose(fid_read);
%         % ----------------------------------
        
        dataall(:,i)=data;
        clear data
    end
    
    datasize = whos('dataall');
    datasize = datasize.bytes;
    
    if datasize < (handles.RAM*10^9)/2
        xy = (dataall' * dataall) / datalen;
    else
        xy = nan(size(dataall,2),'single');
        for i = 1 : size(dataall,1)/chunks
            if i == 1
                xy = dataall(i*chunks-chunks+1:i*chunks,:)' * dataall(i*chunks-chunks+1:i*chunks,:);
            elseif i == size(dataall,1)/chunks
                xy = xy + dataall(i*chunks-chunks+1:end,:)' * dataall(i*chunks-chunks+1:end,:);
            else
                xy = xy + dataall(i*chunks-chunks+1:i*chunks,:)' * dataall(i*chunks-chunks+1:i*chunks,:);
            end
        end
        xy = xy/datalen;
    end
    
    SD = sqrt(diag(xy));
    R = xy./(SD*SD');
    
    COEFF = pcacov(R);
    
    save(['PCAcoeffarraynew' mat2str(j) '.mat'],'xy','SD','R','COEFF')
    
    clear xy R datasize
    
    for i = 1 : ceil(datalen/chunks)
        if i == ceil(datalen/chunks)
            chunksize = size(dataall(i*chunks-chunks+1:end,:),1);
            score = dataall(i*chunks-chunks+1:end,:)./repmat(SD',[chunksize 1]) * COEFF;
            dataFilt = (COEFF(:,max(COEFF) >= thresholdtodel)*score(:,max(COEFF) >= thresholdtodel)').*repmat(SD,[1 chunksize]);
            for k = 1 : length(chanind)
                fwrite(fid_dummy(k + dummyoffset), dataFilt(k,:),'single');
            end
            clear score dataFilt chunksize
        else
            score = dataall(i*chunks-chunks+1:i*chunks,:)./repmat(SD',[chunks 1]) * COEFF;
            dataFilt = (COEFF(:,max(COEFF) >= thresholdtodel)*score(:,max(COEFF) >= thresholdtodel)').*repmat(SD,[1 chunks]);
            for k = 1 : length(chanind)
                fwrite(fid_dummy(k + dummyoffset), dataFilt(k,:),'single');
            end
            clear score dataFilt
        end
    end
    
    dummyoffset = dummyoffset + length(chanind);
    clear dataall COEFF SD chanind
    
    toc
end

for i = 1 : length(chanindall)
    fseek(fid_dummy(i),0, 'bof');
    data = fread(fid_dummy(i),'single=>single');
    fclose(fid_dummy(i));
    delete(['dummy' sprintf('%03d',chanindall(i))]);
    
    data = data';
    save(['datafilt2_ch' sprintf('%03d',chanindall(i))],'data')
    
    clear data
end

clear chanindall datalen
warning on

end



