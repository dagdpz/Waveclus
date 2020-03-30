function Extract_spikes5_tryout(input, varargin)
% 
% files = dir('datafilt2_ch*.mat');
% % files = dir('datafilt_ch*.mat');
% files = {files.name};

%default config
par_input = struct;
parallel = false;

%search for optional inputs
nvar = length(varargin);
for v = 1:2:nvar
    if strcmp(varargin{v},'par')
        if (nvar>=v+1) && isstruct(varargin{v+1})
            par_input = varargin{v+1};
        else
            error('Error in ''par'' optional input.')
        end
    elseif strcmp(varargin{v},'parallel')
        if (nvar>=v+1) && islogical(varargin{v+1})
            parallel = varargin{v+1};
        else
            error('Error in ''parallel'' optional input.')
        end
    end
end

% get a cell of filenames from the input
if isnumeric(input) || any(strcmp(input,'all'))  % cases for numeric or 'all' input
    filenames = {};
    se = supported_wc_extensions();
    dirnames = dir();
    dirnames = {dirnames.name};
    
    for i = 1:length(dirnames)
        fname = dirnames{i};
        [unused, f, ext] = fileparts(fname);
        ext = lower(ext(2:end));
        if any(strcmp(ext,se)) 
            if strcmp(input,'all')
                if strcmp(ext,'mat')
                    warning('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                    continue
                end
                filenames = [filenames {fname}];
            else
                aux = regexp(f, '\d+', 'match');
                if ~isempty(aux) && ismember(str2num(aux{1}),input)
                    if strcmp(ext,'mat')
                        warning('Skipped file ''%s''. The ''.mat'' files should be added by name.\n',fname);
                        continue
                    end
                    filenames = [filenames {fname}];   
                end
            end
        end
    end
    
elseif ischar(input) && length(input) > 4  
    if  strcmp (input(end-3:end),'.txt')   % case for .txt input
        filenames =  textread(input,'%s');
    else
        filenames = {input};               % case for cell input
    end
    
elseif iscellstr(input)   %case for cell input
    filenames = input;
else
    ME = MException('MyComponent:noValidInput', 'Invalid input arguments');
    throw(ME)
end


    par = set_parameters();
    par.filename = filenames{1};
    par.reset_results = true;  %if true,  don't load times_ or _spikes files
    par.cont_segment = true;  %false doesn't save the segment of the continuous data in the spikes file
    try 
        data_handler = readInData(par);
    catch MExc
        warning(MExc.message);
        return
    end
    
    par = data_handler.par;
    par = update_parameters(par,par_input,'detect');
    data_handler.par = par;

sr = par.sr;
w_pre = par.w_pre;
w_post = par.w_post;
ref = par.ref;
detect = par.detection;
% stdmin = par.stdmin;
% stdmax = par.stdmax;
% fmin_detect = par.detect_fmin;
% fmax_detect = par.detect_fmax;
% fmin_sort = par.sort_fmin;
% fmax_sort = par.sort_fmax;

for k= 1:length(filenames)
    tic
    load(filenames{k},'data')
    
%    data = double(data);
%     if ~strcmp(handles.dtypewrite,'single')
%         data=single(data);
%     end
    
    if ~(par.transform_factor == 1)
        data = data * par.transform_factor;                       %%in microvolts
    end
    
    if size(data,2) == 1
        data = data';
    end
    
    numsamples = length(data);
    
    thr = par.stdmin*median(abs(data))/0.6745;
    thrmax = par.stdmax*median(abs(data))/0.6745;
    
    % LOCATE SPIKE TIMES
    switch par.detection
        case 'pos', ups = find(data(w_pre+2:end-w_post-2) > thr) +w_pre+1;%indices of values above threshold
        case 'neg', ups = find(data(w_pre+2:end-w_post-2) < -thr) +w_pre+1;%indices of values below threshold
        case 'both', ups = find(abs(data(w_pre+2:end-w_post-2)) > thr) +w_pre+1;   
    end
    if isempty(ups), %no spikes detected
        spikes=[]; index=[];
        %    save(spikes_fname,'spikes','index');
    else
        
        %find beginnings and endings of intervals with values above threshold
        begs=ups([true diff(ups)>1]);%indices in data array
        ends=ups([diff(ups)>1 true]);
        %find two events which are closer than 1.4ms, remove smallest from two
        %closest
        bad=[];
        for i=1:length(begs)-1,
            if begs(i+1)-ends(i)<=ref*sr 
                bad=[bad i]; 
            end
        end
        bad2=[];
        for i=1:length(begs)-2,
            if begs(i+2)-ends(i)<=ref*sr 
                bad2=[bad2 i]; 
            end
        end

        nspk=length(begs);
        index=nan(nspk,1);
        np=w_post+w_pre;
        spikes=nan(nspk,np+4);
        m=nan(nspk,1);
        
        for i=1:nspk,
            [m(i),ind] = max(abs(data(begs(i):ends(i))));
            index(i) = begs(i) + ind -1;
            spikes(i,:)= data(index(i)-w_pre+1-2:index(i)+w_post+2);
        end
        clear data
        badind=[];
        if ~isempty(bad)
            for i=bad,
                if m(i)>m(i+1)
                    badind=[badind i+1];
                else
                    badind=[badind i];
                end
            end
        end
        if ~isempty(bad2)
            for i=bad2,
                if m(i)>m(i+2)
                    badind=[badind i+2];
                else
                    badind=[badind i];
                end
            end
        end
        badind = unique(sort(badind));
        ind=setdiff(1:nspk,badind);
        ind = ind(m(ind) <= thrmax);
        spikes=spikes(ind,:);
        index=index(ind);
        nspk=length(ind);
        %INTERPOLATION
        if par.interpolation=='y',
            %Does interpolation
            ints=1/int_factor:1/int_factor:np+4;
            intspikes = spline(1:np+4,spikes,ints);
            spikes1=zeros(nspk,np*int_factor);
            
            switch par.detection
                case 'pos'
                    [~,iaux]=max(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
                case 'neg'
                    [~,iaux]=min(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
                case 'both'
                    [~,iaux]=max(abs(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor)),[],2);
            end
            iaux = iaux + (w_pre+1)*int_factor -1;
            for i=1:nspk
                spikes1(i,w_pre*int_factor:-1:1) = intspikes(i,iaux(i):-1:iaux(i)-(w_pre*int_factor-1));
                spikes1(i,w_pre*int_factor+1:np*int_factor) = intspikes(i,iaux(i)+1:1:iaux(i)+w_post*int_factor);
            end
            spikes=spikes1;
            clear spikes1
            clear intspikes
        else
            spikes(:,[1 2 end-1 end])=[];       %eliminates borders that were introduced for interpolation
        end
        
        index=index/sr*1000;
        
%         chanceevent = numsamples * chancelevel/2;    
        
        switch par.detection
            case 'pos'
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr')
%                 end
            case 'neg'
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_onethr'],'spikes','index','thr')
%                 end
            case 'both'
                spikesign = sign(spikes(:,w_pre*int_factor));
                spikesall = spikes;
                indexall = index;
                spikes = spikesall(spikesign == -1,:);
                index = indexall(spikesign == -1);
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_negthr'],'spikes','index','thr')
%                 end
                spikes = spikesall(spikesign == 1,:);
                index = indexall(spikesign == 1);
%                 if length(index) > chanceevent
                    save(['dataspikes_ch' files{k}(strfind(files{k},'ch')+2:end-4) '_posthr'],'spikes','index','thr')
%                 end
        end
        
       
    end
    clearvars -except handles sr w_pre w_post int_factor files k chancelevel
    toc
end


