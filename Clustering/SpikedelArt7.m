function SpikedelArt7(badchan,handles)

chunks = 1000000;
thresholdtodel = 0.36;

warning off

for j = 1 : handles.numArray
    tic
    k = 1;
    fid_read=[];
    for cc = 1 : handles.numchan
        channel=['datahpch' sprintf('%03d',(cc+j*handles.numchan-handles.numchan))];
        if sum(badchan == str2double(channel(end-2:end)))
        else
            readName = ['datahpch' sprintf('%03d', (cc+j*handles.numchan-handles.numchan))];
            fid_read(k) = fopen(readName);
            k = k+1;
        end
        clear readName channel
    end
    clear k
    
    data = fread(fid_read(1),'int16=>int16');
    dataall=zeros(length(data),length(fid_read),'int16');
    dataall(:,1)=data;
    for i = 2 : length(fid_read)
        data = fread(fid_read(i),'int16=>int16');
        dataall(:,i)=data;
        clear data
    end
    
    for cc = 1:length(fid_read)
        fclose(fid_read(cc));
    end
    
    [m,n] = size(dataall);
    
    for cc = 1:n
        fid_write(cc) = fopen(['dummy' sprintf('%03d', cc )],'w');
    end
    
    dataall = single(dataall);

  % direct calculation of the covariance matrix
%     xy = (dataall' * dataall) / m;

  % indirect but memory saving version for calculation of the covariance matrix  
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
    xy = xy/m;
    
    SD = sqrt(diag(xy));
    R = xy./(SD*SD');
    
    COEFF = pcacov(R);
    
    save(['PCAcoeffarray' mat2str(j) '.mat'],'xy','SD','R','COEFF')
    
    clear xy R 
    
    for i = 1 : ceil(m/chunks)
        if i == ceil(m/chunks)
            chunksize = size(dataall(i*chunks-chunks+1:end,:),1);
            score = COEFF\(dataall(i*chunks-chunks+1:end,:)./repmat(SD',[chunksize 1]))';
            dataFilt = (COEFF(:,max(COEFF) >= thresholdtodel)*score(max(COEFF) >= thresholdtodel,:)).*repmat(SD,[1 chunksize]);
            for cc = 1:n
                fwrite(fid_write(cc), dataFilt(cc,  :) , 'int16');
            end
            clear score dataFilt chunksize
        else
            %         score = inv(COEFF)*dataall(i*chunks-chunks+1:i*chunks,:)';
            score = COEFF\(dataall(i*chunks-chunks+1:i*chunks,:)./repmat(SD',[chunks 1]))';
            dataFilt = (COEFF(:,max(COEFF) >= thresholdtodel)*score(max(COEFF) >= thresholdtodel,:)).*repmat(SD,[1 chunks]);
            for cc = 1:n
                fwrite(fid_write(cc), dataFilt(cc,  :) , 'int16');
            end
            clear score dataFilt
        end
    end
    
    clear dataall COEFF SD
    
    for cc = 1:n
        fclose(fid_write(cc));
    end
    
    counter = 1;
    for i = 1 : handles.numchan
        channel=['datahpch' sprintf('%03d',(i+j*handles.numchan-handles.numchan))];
        if sum(badchan == str2double(channel(end-2:end)))
        else
            readName = ['dummy' sprintf('%03d', counter)];
            fid_read = fopen(readName);
            data = fread(fid_read,'int16=>int16');
            fclose(fid_read);
            
            save(['datahighpassch' sprintf('%03d',i+j*handles.numchan-handles.numchan)],'data')
            
            counter = counter+1;
            clear data readName
            
        end
        clear channel
    end
    
    for cc = 1:n
        delete(['dummy' sprintf('%03d',cc)]);
    end
    
    clear m n channel counter
    
    toc
    
end

warning on


