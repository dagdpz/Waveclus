function [inspk,feature_names,inputs] = wave_features6(spikes,handles)
%Calculates the spike features

exclusioncrit = handles.par.exclusioncrit; % thr; number 
exclusionthr = handles.par.exclusionthr;
maxinputs = handles.par.maxinputs;
int_factor = handles.par.int_factor;
w_pre = handles.par.w_pre;
w_post = handles.par.w_post;
% sr = handles.par.sr;

shift = 2;
offstart = round(w_pre/2) + shift;
offend  = round(w_post/2) + shift + w_pre;

ind1 = [fliplr(w_pre : -1: offstart+1),w_pre+1 : offend] * int_factor;

% downfac = 2;
% ind2 = [fliplr(w_pre : -downfac : offstart+1),w_pre+1 : downfac : offend] * int_factor;

% CALCULATES FEATURES

%% PCA
    
[~,S,latent] = princomp(spikes(:,ind1));

switch exclusioncrit
    
    case 'thr'
        latent = cumsum(latent./sum(latent));
        
        inputs = find(latent > exclusionthr,1,'first');
        
    case 'number'
        inputs = maxinputs;
        
end
        
%features
feature_names = cell(1,inputs);
for i=1:inputs
    feature_names{i}=sprintf('PCA,%d',i);
end
        
inspk = S(:,1:inputs);

fprintf('Selected freatures: ');
fprintf('%s ',feature_names{:});
fprintf('\n')

% %CREATES INPUT MATRIX FOR SPC
% inspk=zeros(nspk,inputs);
% for i=1:nspk
%     for j=1:inputs
%         inspk(i,j)=cc(i,coeff(j));
%     end
% end

