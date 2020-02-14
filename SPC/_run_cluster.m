function [clu, tree] = run_cluster(handles);
dim=handles.par.inputs;
fname=handles.par.fname;
%set(handles.file_name,'string','Running SPC ...');

% DELETE PREVIOUS FILES
save([fname '.dg_01.lab'],'dim','-ASCII');          delete([fname '.dg_01.lab']);
save([fname '.dg_01'],'dim','-ASCII');              delete([fname '.dg_01']);

dat=load(fname);
%dat=load('data');
n=length(dat);
fid=fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname);
fprintf(fid,'Dimentions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(handles.par.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(handles.par.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(handles.par.tempstep));
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'SWCycles: %s\n',num2str(handles.par.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(handles.par.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
fclose(fid);

switch handles.par.system
    case 'windows'
        if exist([pwd '\cluster.exe'])==0
            directory = which('cluster.exe');
            copyfile(directory,pwd);
        end
        dos(sprintf('cluster.exe %s.run',fname));
%        !cluster data.run
    case 'linux'
        if exist([pwd '/cluster_linux.exe'])==0
            directory = which('cluster_linux.exe');
            copyfile(directory,pwd);
        end
        unix(sprintf('./cluster_linux %s.run',fname));
end
        
clu=load([fname '.dg_01.lab']);
tree=load([fname '.dg_01']); 
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname);