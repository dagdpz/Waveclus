function [inspk,feature_names] = wave_features_new(spikes,handles)
%Calculates the spike features

scales = handles.par.scales;
feature = handles.par.features;
inputs = handles.par.inputs;
nspk=size(spikes,1);
len = size(spikes,2);
%set(handles.file_name,'string','Calculating spike features ...');

% CALCULATES FEATURES
switch feature
    case 'wav'
        cc=zeros(nspk,len);
        for i=1:nspk                                % Wavelet decomposition
            [c,l]=wavedec(spikes(i,:),scales,handles.par.wavelet);
            cc(i,1:len)=c(1:len);
        end
        %wavelet features
        fn=cell(1,len);
        j=1;
        k=scales+1;
        ll=cumsum(l);
        for i=ll(1:end-1),
            j1=1;
            while j<=i, 
                if i==ll(1), fn{j}=sprintf('A%d,%d',k-1,j1);
                else fn{j}=sprintf('D%d,%d',k,j1); end
                j=j+1; j1=j1+1; 
            end
            k=k-1;
        end
    case 'pca'
        [C,S] = princomp(spikes);
        cc = S;
        inputs = 3;
        %pca features
        fn=cell(1,3);
        for i=1:3, fn{i}=sprintf('PCA,%d',i); end
   case 'wavpca',
      %wavelet
       cc=zeros(nspk,len);
            for i=1:nspk                                % Wavelet decomposition
                [c,l]=wavedec(spikes(i,1:len/2),scales,handles.par.wavelet);
                cc(i,1:len/2)=c(1:len/2);             %def.cc(i,1:ls)=c(1:ls);
                [c,l]=wavedec(spikes(i,len/2+1:end),scales,handles.par.wavelet);
                cc(i,len/2+1:len)=c(1:len/2);
            end
        %pca
      [C,S] = princomp(spikes(:,1:len/2));
      cc(1:nspk,len+1:len+3) = S(:,1:3);
      [C,S] = princomp(spikes(:,len/2+1:len));
      cc(1:nspk,len+4:len+6) = S(:,1:3);
      
      %assign names of features
      fn=cell(1,len+6);
      %wavelet features
      j=1;
      k=scales+1;
      ll=cumsum(l);
      for i=ll(1:end-1),
         j1=1;
         while j<=i,
            if i==ll(1), fn{j}=sprintf('A%d,%d',k-1,j1);
                fn{j+len/2}=sprintf('AX%d,%d',k-1,j1);
            else fn{j}=sprintf('D%d,%d',k,j1);
                fn{j+len/2}=sprintf('DX%d,%d',k,j1);
            end
            j=j+1; j1=j1+1;
         end
         k=k-1;
      end
      %pca features
      for i=1:3, fn{len+i}=sprintf('PCA,%d',i); 
        fn{len+3+i}=sprintf('PCAX,%d',i); 
      end

end
warning('off','stats:lillietest:OutOfRangeP');
nf=size(cc,2);
sd=NaN*ones(1,nf);
for i=1:nf,                            % KS test for coefficient selection
   thr_dist = std(cc(:,i)) * 3;
   thr_dist_min = mean(cc(:,i)) - thr_dist;
   thr_dist_max = mean(cc(:,i)) + thr_dist;
         aux = cc((cc(:,i)>thr_dist_min)&(cc(:,i)<thr_dist_max),i);
%          ind=find((cc(:,i)<=thr_dist_min)|(cc(:,i)>=thr_dist_max));
%          cc(ind,i)=NaN;
%          aux = cc(:,i);

%   warning off
   [hh,pp,lstat]=lillietest(aux);
%    warning on
   sd(i)=lstat;
end
% fprintf('Lstat: ');
% fprintf('%1.2f ',sd);
% fprintf('\n');
warning('on','stats:lillietest:OutOfRangeP');
[max ind]=sort(-sd);
feature_names={fn{ind(1:inputs)}};
inspk=cc(:,ind(1:inputs));
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

