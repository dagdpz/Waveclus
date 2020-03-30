function [clu, tree] = wc_run_cluster(handles)
dim=handles.WC.inputs;
fname=handles.fname;

% DELETE PREVIOUS FILES, if exist
if exist([fname '.dg_01.lab'],'file'), delete([fname '.dg_01.lab']); end
if exist([fname '.dg_01'],'file'), delete([fname '.dg_01']); end

dat=load(fname);
%dat=load('data');
n=length(dat);
fid=fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname);
fprintf(fid,'Dimentions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(handles.WC.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(handles.WC.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(handles.WC.tempstep));
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'SWCycles: %s\n',num2str(handles.WC.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(handles.WC.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
fclose(fid);


if ispc, sys = 'windows'; 
   cluster_prog='cluster.exe';
elseif ismac, sys = 'MACI64'; 
   cluster_prog='cluster_maci.exe';
else
   if isunix, sys = 'linux';
      cluster_prog='cluster_linux';
   else error('run_cluster:UnknownSystem','Unknown system');
   end
end
loc = which(cluster_prog);
if isempty(loc), error('run_cluster:clusterprogNotFound','SPC programm is not found (cluster.exe or cluster_linux.exe)'); end

switch sys
    case 'windows'
        dos(sprintf('%s %s.run',loc,fname));
    case 'linux'
        unix(sprintf('%s %s.run',loc,fname));
    case 'MACI64'
        unix(sprintf('%s %s.run',loc,fname));
end
        
clu=load([fname '.dg_01.lab']);
tree=load([fname '.dg_01']); 
%delete temporaray files
delete([fname '.dg_01.lab']);
delete([fname '.dg_01']);
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname);