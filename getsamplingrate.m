function [sr]=getsamplingrate(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 IMPORT OF SPIKE EVENTS                
%
% DESCRIPTION:
% This routine returns the sampling rate of a recording system. Function
% retursn -1 if sampling rate was not found.
% 
% HELPFUL INFORMATION:  -How to import TDT Tank into Matlab.pdf
%                       -Tank format.pdf
%
% SYNTAX:    [sr] = getsamplingrate(system,'TDT',path,'path');
%
%        system    ... char(identifier of system;'TDT' Tucker Davis Techn.
%                                                'BR'  Blackrock Microsys.
%        path     ... path of recording
%        sr        ... sampling rate in Hz [1/s]
%                       
%
% EXAMPLE:   
%          sr = getsamplingrate('path','C:\TANK1\Block-16','system','TDT')
%
% AUTHOR: ©Stefan Schaffelhofer, German Primate Center              AUG11 %
%
%          DO NOT USE THIS CODE WITHOUT AUTHOR'S PERMISSION! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
% Error check:

if nargin<4
    error('Not enough input arguements.');
end
if nargin>4
    error('Too many input arguements.');
end

pathok   = 0;
systemok = 0;


% Check if all necessary input parameters are there
for ii=1:numel(varargin)
    if strcmpi(varargin{ii},'path')
        pathok=1;
        pathid=ii+1;
    elseif strcmpi(varargin{ii},'system')
        systemok=1;
        systemid=ii+1;
    end
    
end

if ~(pathok && systemok)
    error('Not enought input arguements. Path and import-system must be defined.');
end

% Check if system is known
if isa(varargin{systemid},'char')
    if strcmpi(varargin{systemid},'TDT')                                   % Tucker Davis Technologies
        system='TDT';
    elseif strcmpi(varargin{systemid},'BL')                                % Blackrock
        system='BL';
    else
        error('Unknown system.')
    end
else
    error(['Wrong option of ' varargin{systemid-1} '- parameter.']); 
end

% Check if folder/file, or tank exists
if isa(varargin{pathid},'char')
    switch system
        case 'TDT'                                                         % Check if folders/files exist
                if exist(varargin{pathid},'dir')
                   thispath=varargin{pathid};
                   tsqfiles=dir([thispath '/*.tsq']);
                   tevfiles=dir([thispath '/*.tev']);
                   if numel(tsqfiles)==1 && numel(tevfiles)==1
                       tsq_path=[thispath '/' tsqfiles.name];
                   else
                       error('Two or more *.tsq or *.tev files.');
                   end
                else
                    error([system '-tank not found.']);
                end
        case 'BL'
    end

end

% IMPORT data
switch system
    case 'TDT'                                                             % Import for Tucker Davis Technologies
    % open the files
    tsq = fopen(tsq_path,'r');                                             % ... and tell the position (=number of bytes from beginning). Dividing the number of bytes of the file by the number of bytes for a header blocks gives back the number of headers.
    fseek(tsq, 0, 'bof');                                                  % go back to the start of the file
    disp('Load TSQ file ...');
    namesearcheNeu = 256.^(0:3)*double('Wave')';
    %namesearcheNeu = 256.^(0:3)*double('eNeu')';
    %namesearcheNeu = 256.^(0:3)*double('pNeu')';
    k=0;
    sr=-1;
    while k==0
        TSQ=fread(tsq,40,'uint8=>uint8');
        if ~isempty(TSQ)
            name=typecast(TSQ(9:12),'uint32');
            if name==namesearcheNeu
                sr = double(typecast(TSQ(37:40,1),'single'));
                k=1;
            end
        else
            k=1;                                                           % no spikes found, go out of loop, retrun -1;
        end
    end
    
    fclose(tsq);                                                           % close files

    
    case 'BL'                                                              % Import for Blackrock
end

t=toc;
disp(['Spikes imported in ' num2str(t) ' sec.']);
end
