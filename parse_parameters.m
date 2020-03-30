function p=parse_parameters(v,var_names,defaults)

% var_names={'ints','bin','tmax','title','color','norm'};
% defaults={[],1,100,'Autocorrelogram','k',0};
for i=1:length(var_names),
   p.(var_names{i})=defaults{i};
end
if floor(length(v)/2)~=length(v)/2, v(end)=[];end
for i=1:2:length(v),
   s=v{i};
   if ischar(s), 
      t=strmatch(lower(s),var_names,'exact');
      if ~isempty(t),
         p.(var_names{t})=v{i+1};
      end
   else
      fprintf('parse_parameters:Something is wrong with optional arguments, just ignoring them\n');
   end
end
end