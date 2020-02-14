function SpikefilterChan8(handles)

stop_f = handles.linenoisefrequ;
sr = handles.par.sr;
hpcutoff = handles.hpcutoff;
lpcutoff = handles.lpcutoff;
transform_factor = handles.par.transform_factor;
rawname = handles.rawname;
dtyperead = handles.dtyperead;
dtypewrite = handles.dtypewrite;

if strcmp(handles.hp,'med')
    n = floor(sr/hpcutoff);
    if ~mod(n,2)
        n = n+1;
    end
end   


if handles.blockfile
    files = dir('dataraw_ch*.mat');
    files = {files.name};
else
    files = dir(rawname);
    files = {files.name};
end

warning off

for i = 1 : length(files)
    tic
    
    switch handles.sys
        case 'TD'
            
            filename =  files{i};
            channel_file_handle = fopen(filename, 'r');                       %'r' - read the indicated channel file
            fseek(channel_file_handle, 0, 'eof');                                      %sets the file pos indicator at the end of the channel file
            store = ftell(channel_file_handle);                                        %gives the total number of bytes in the channel file
            
            time_window(1) = 0;
            time_window(2) = (store/2)/sr;                              %in sec
            
            no_samples = round(sr*(time_window(2) - time_window(1)));              %number of samples during the indicated time period (sampling rate = samples per ms)
            start_ind = sr*(time_window(1));                                %where to start to read out the data
            
            fseek(channel_file_handle, start_ind*2, 'bof');                            %sets the file pos indicator at the byte where the read out is supposed to start
            data = fread(channel_file_handle, no_samples, [dtyperead '=>' dtyperead]);             %reads out the data, int16 - class of input and output values
            fclose(channel_file_handle);
            
            data = double(data)/transform_factor;
            
            if handles.iniartremovel
                data(1:40) = mean(data(41:80));
            end
            
        case 'BR'
            load(files{i});
            data = double(data);
            
        case 'RHD2000'
            load(files{i});
            data = double(data);
    end
    
    if size(data,2) == 1
        data = data';
    end
    
    namefile = files{i}
    
    if handles.linenoisecancelation
        
        [b,a]=ellip(2,0.1,40,[stop_f-1 stop_f+1]*2/sr,'stop');
        data=filtfilt(b,a,data);
        [b,a]=ellip(2,0.1,40,[2*stop_f-1 2*stop_f+1]*2/sr, 'stop');
        data=filtfilt(b,a,data);
        [b,a]=ellip(2,0.1,40,[3*stop_f-1 3*stop_f+1]*2/sr, 'stop');
        data=filtfilt(b,a,data);
       
    end
    
    if strcmp(handles.hp,'med')
        xx = median_filter(data,n);
        xxend = medfilt1(data(end-n:end),n);
        xx = [xx(ceil(n/2):end) xxend(ceil(n/2)+2:end)];
        data = data-xx;
        clear xxend xx;
    else
        [b,a] = butter(2,hpcutoff*2/sr,'high');
        data=filtfilt(b,a,data);
    end
    
        [b,a] = butter(2,lpcutoff*2/sr,'low');
        data=filtfilt(b,a,data);
        
        data = eval([dtypewrite '(data)']);
        save(['datafilt_ch' sprintf('%03d',str2double(namefile(strfind(namefile,'ch')+2:end-4)))],'data');
    
    clear namefile data
    toc
end
clear files

warning on



