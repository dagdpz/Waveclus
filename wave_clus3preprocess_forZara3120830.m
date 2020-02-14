clear all
dbstop if error

count = 1;

% RECORDING SYSTEM
% handles.sys = 'TD';
handles.sys = 'BR';
% handles.sys = 'SpikeSim';

% SAMPLING RATE
% handles.par.sr = 24414.0625;
handles.par.sr = 30000;
% handles.par.sr = 25000;

% ARRAY CONFIGURATION
handles.numArray = 4;    
handles.numchan = 32;         %Channels per Array


% LINE NOISE
handles.linenoisecancelation = 0;      % 1 for yes; 0 for no
% handles.linenoisefrequ = 100;
handles.linenoisefrequ = 50;

% FILTER OPTIONS
handles.hp = 'med';                    % med = medianfiltersubtraction; int = interpolationsubtraction; but = butterworth
handles.hpcutoff = 333;                % in Hz
handles.lpcutoff = 5000;               % in Hz
handles.par.transform_factor = 0.25;      % microVolts per bit for higher accuracy when saved as int16 after filtering; Default for BR
handles.iniartremovel = 0;
handles.drinkingartremoval = 0;

% DETECTION
handles.par.w_pre = 20;
handles.par.w_post = 44;
handles.par.int_factor = 2;
handles.par.interpolation ='y';
handles.par.stdmin = 5;
handles.par.stdmax = 100; 
handles.threshold ='both';   


% Zara 
badchan = {[]};



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
            SpikefilterChan7(handles);
            SpikedelArt7(badchan{count},handles); 
%             for cc = 1 : handles.numchan*handles.numArray
%                 delete(['datahpch' sprintf('%03d',cc)]);
%             end
            txt_files = exist('files.txt','file');
            if(txt_files ~=0 )
                delete('files.txt');
            end
            if handles.drinkingartremoval
                Spikecuttrials2
                matlist = dir('datacut*.mat');
            else
                matlist = dir('datahighpassch*.mat');
            end
            for i=1:length(matlist);
                name=matlist(i).name(1:end-4);
                dlmwrite('files.txt',name,'-append','delimiter','')
            end
            Extract_spikes2(handles);
            Do_clustering2_redo(handles);
            clear functions
            count = count + 1;
        end
        disp('This is the end ...')
    end
end