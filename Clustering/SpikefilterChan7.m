function SpikefilterChan7(handles)

stop_f = handles.linenoisefrequ;
sr = handles.par.sr;
hpcutoff = handles.hpcutoff;
lpcutoff = handles.lpcutoff;
transform_factor = handles.par.transform_factor;

if strcmp(handles.hp,'med')
    n = floor(sr/hpcutoff);
    if ~mod(n,2)
        n = n+1;
    end
end   

mex /Users/benjaminwellner/MatTool/functions/median_filter.c

switch handles.sys
    case 'TD'
        files = dir('Moe_Recording*.sev');
        filesB = {files.name};
    case 'BR'
        files = dir('*_c_ch*.mat');
        filesB = {files.name};
    case 'SpikeSim'
        filesB = {'datafile_ch001.mat';'datafile_ch002.mat';'datafile_ch003.mat'};
end


warning off

for i = 1 : length(filesB)
    tic
    
    switch handles.sys
        case 'TD'
            
            filename =  filesB{i};
            channel_file_handle = fopen(filename, 'r');                       %'r' - read the indicated channel file
            fseek(channel_file_handle, 0, 'eof');                                      %sets the file pos indicator at the end of the channel file
            store = ftell(channel_file_handle);                                        %gives the total number of bytes in the channel file
            
            time_window(1) = 0;
            time_window(2) = (store/2)/sr;                              %in sec
            
            no_samples = round(sr*(time_window(2) - time_window(1)));              %number of samples during the indicated time period (sampling rate = samples per ms)
            start_ind = sr*(time_window(1));                                %where to start to read out the data
            
            fseek(channel_file_handle, start_ind*2, 'bof');                            %sets the file pos indicator at the byte where the read out is supposed to start
            data = fread(channel_file_handle, no_samples, 'int16=>int16');             %reads out the data, int16 - class of input and output values
            if size(data,2) == 1
                data = data';
            end
            fclose(channel_file_handle);
            
            x = double(data)/transform_factor;
            clear data
            
            if handles.iniartremovel
                x(1:40) = mean(x(41:80));
            end
            
        case 'BR'
            load(filesB{i});
            
            x = double(data);
            clear data
            
        case 'SpikeSim'
            files = dir('ST2013*.mat');
            load(files.name,'signals');
            x = signals(:,i)'/transform_factor;
            clear signals
    end
    
    namefile = filesB{i}
    
    fid_write = fopen(['datahpch' sprintf('%03d',str2double(namefile(strfind(namefile,'ch')+2:end-4)))],'w');
    
    if handles.linenoisecancelation
        
        [b,a]=ellip(2,0.1,40,[stop_f-1 stop_f+1]*2/sr,'stop');
        x=filtfilt(b,a,x);
        [b,a]=ellip(2,0.1,40,[2*stop_f-1 2*stop_f+1]*2/sr, 'stop');
        x=filtfilt(b,a,x);
        [b,a]=ellip(2,0.1,40,[3*stop_f-1 3*stop_f+1]*2/sr, 'stop');
        x=filtfilt(b,a,x);
        
%         [b,a]=ellip(2,0.1,40,[4*stop_f-2 4*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         [b,a]=ellip(2,0.1,40,[5*stop_f-2 5*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         [b,a]=ellip(2,0.1,40,[7*stop_f-2 7*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         
%         [b,a]=ellip(2,0.1,40,[9*stop_f-2 9*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         
%         [b,a]=ellip(2,0.1,40,[11*stop_f-2 11*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         [b,a]=ellip(2,0.1,40,[15*stop_f-2 15*stop_f+2]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         
%         [b,a]=ellip(2,0.1,40,[5130 5170]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
%         
%         [b,a]=ellip(2,0.1,40,[6600 6630]*2/sr, 'stop');
%         x=filtfilt(b,a,x);
        
    end
    
    if strcmp(handles.hp,'med')
        xx = median_filter(x,n);
        xxend = medfilt1(x(end-n:end),n);
        xx = [xx(ceil(n/2):end) xxend(ceil(n/2)+2:end)];
        x = x-xx;

%         [b,a] = butter(4,hpcutoff*2/sr,'high');
%         x=filtfilt(b,a,x);

        clear xxend xx;
    elseif strcomp(handles.hp,'int')
        % not done jet
    else
        [b,a] = butter(2,hpcutoff*2/sr,'high');
        x=filtfilt(b,a,x);
    end
    
        [b,a] = butter(2,lpcutoff*2/sr,'low');
        x=filtfilt(b,a,x);
    
%         [b,a] = butter(8,lpcutoff*2/sr,'low');
%         x=filtfilt(b,a,x);
        
        fwrite(fid_write, x , 'int16');
        clear xx

    fclose(fid_write);
    
    clear namefile
    toc
end
clear files filesB

warning on



