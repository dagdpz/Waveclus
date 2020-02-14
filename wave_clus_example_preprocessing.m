clear all
dbstop if error

%mex(which('median_filter.c'))
% mex /Users/benjaminwellner/MatTool/functions/median_filter.c

count = 1;

% SYSTEM MEMORY
handles.RAM = 64;                                                           % in GB

% RECORDING SYSTEM
handles.sys = 'TD';
% handles.sys = 'BR';
% handles.sys = 'RHD2000';

% RAW DATAFILES NAME
% handles.rawname = '*_c_ch*.mat';
handles.rawname = '*_Recording*.sev';
% handles.rawname = 'amplifier*.bin';

handles.blockfile = 0;


% SAMPLING RATE
handles.par.sr = 24414.0625;

handles.dtyperead = 'int16';                                                  % Default for BR, TD

if strcmp(handles.dtyperead,'float32')
    handles.dtype = 'single';
end
if strfind(handles.dtyperead,'uint');
    handles.dtypewrite = handles.dtyperead(2:end);
else
    handles.dtypewrite = handles.dtyperead;
end

% ARRAY CONFIGURATION
handles.numArray = 6;    
handles.numchan = 32;                                                       % Channels per Array
handles.arraynoisecancelation = 1;

% LINE NOISE
handles.linenoisecancelation = 0;                                           % 1 for yes; 0 for no
% handles.linenoisefrequ = 100;
handles.linenoisefrequ = 50;

% FILTER OPTIONS
handles.hp = 'med';                                                         % med = medianfiltersubtraction; int = interpolationsubtraction; but = butterworth
handles.hpcutoff = 333;                                                     % in Hz
handles.lpcutoff = 5000;                                                    % in Hz
% handles.par.transform_factor = 0.25;                                      % microVolts per bit for higher accuracy when saved as int16 after filtering; Default for BR
handles.par.transform_factor = 0.25;                                           % Default for RHD2000
handles.iniartremovel = 1;
handles.drinkingartremoval = 0;

% DETECTION
handles.par.w_pre = 20;
handles.par.w_post = 44;
handles.par.int_factor = 2;
handles.par.interpolation ='y';
handles.par.stdmin = 5;
handles.par.stdmax = 100; 
handles.threshold ='both';   


% Moe 
badchan = {
    [16 26 61 67 71 113 173]
    [16 26 61 67 71 113 173]
    [16 26 61 67 71 113 173]
    };
    
tankname='20170317';
tank = ['T:\Data\TDTtanks\Linus_phys\' tankname '\'];                                                   % path of recordings
recordingnames={'Block-1'};                                            % recording names


%% need to switch to the respective folder...
current_path=pwd;
tankname=num2str(tankname);
for ii =1:length(recordingnames)
    parse_data_tdt([tank recordingnames{ii}])
    
    folder=['WC_' recordingnames{ii} '\'];
    block=recordingnames{ii}(strfind(recordingnames{ii},'-')+1:end);
    recname=['blocks_' block '-01'];
    disp(['Processing: ' recname]);
    cd([tank folder])
    
    
            SpikefilterChan8(handles);
    
    
            if handles.arraynoisecancelation
            SpikedelArt8(badchan{count},handles); 
            end

            Extract_spikes3(handles);
            
            Do_clustering4_redo(handles);
            
            count = count + 1;
            
            
%     cd([tank folder])
%     Get_spikes('all')
%     Do_clustering('all')
%     cd(current_path)
    
    SPK = WC2SPK([tank folder]);                   % convert wave clus files to spike object and neuronobj
    SPK2NEX(SPK,[tank tankname '_' recname '.nex']);
   
    save([tank tankname '_' recname '_SPK.mat'],'SPK');                % save spike object
end
    
load_dir_start = uigetdir(pwd, 'Choose your directory...');
if load_dir_start ~= 0
    load_folders = dir(load_dir_start);
    load_folders = {load_folders([load_folders.isdir]).name};
    load_folders = setxor(load_folders,{'.','..'});
    [bert, ok] = listdlg('ListString',load_folders);
  
    if ok
        for lauf = bert            
            load_dir = fullfile(load_dir_start,load_folders{lauf});
            cd(load_dir)
            
            if handles.blockfile
            Split_Chan(handles);
            end
           
            SpikefilterChan8(handles);
            
            if handles.arraynoisecancelation
            SpikedelArt8(badchan{count},handles); 
            end

            Extract_spikes3(handles);
            
            Do_clustering4_redo(handles);
            
            count = count + 1;
        end
        disp('Done')
    end
end