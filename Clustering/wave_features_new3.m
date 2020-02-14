function [inspk,feature_names,spikesadd] = wave_features_new3(spikes,handles)
%Calculates the spike features

scales = handles.par.scales;
feature = handles.par.features;
inputs = handles.par.inputs;
nspk=size(spikes,1);
len = size(spikes,2);
pcaoffset = 10;
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
        [C,S] = princomp(cat(2,spikes(:,1+pcaoffset:len/2-pcaoffset),spikes(:,len/2+1+pcaoffset:len-pcaoffset)));
        cc(1:nspk,len+1:len+3) = S(:,1:3);
        clear C S
        %integral
        cc(1:nspk,len+4:len+19) = cumsum(spikes(:,15:30)')';
        
        %peaktopeak
        [b,a] = butter(2,3000/15000,'low');
        first = filtfilt(b,a,spikes(:,20:64)');
        seco = filtfilt(b,a,spikes(:,84:end)');
        firstB = diff(first(5:end,:))<=0;
        firstB = logical(cat(1,firstB,ones(1,size(first,2))));
        secoB = diff(seco(5:end,:))<=0;
        secoB = logical(cat(1,secoB,ones(1,size(first,2))));
        peakindfirst = zeros(size(first,2),1);
        firstamp = zeros(size(first,2),1);
        peakindseco = zeros(size(seco,2),1);
        secoamp = zeros(size(seco,2),1);
        for j = 1 : size(first,2)
            dummy = 5:size(first,1);
            peakindfirst(j) = min(dummy(firstB(:,j)));
            clear dummy
            firstamp(j) = first(peakindfirst(j),j)-spikes(j,20)';
            dummy = 5:size(seco,1);
            peakindseco(j) = min(dummy(secoB(:,j)));
            clear dummy
            secoamp(j) = seco(peakindseco(j),j)-spikes(j,84)';
        end
        
        clear seco secoB first firstB b a 
        cc = cat(2,cc,peakindfirst,firstamp,peakindseco,secoamp);
        clear peakindfirst firstamp peakindseco secoamp
        %assign names of features
        %       fn=cell(1,len+6);
        spikesadd = cc(:,len+1:end);
        
        fn=cell(1,len+23);
        
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
        end
        %integral features
        for i=1:16, fn{len+3+i}=sprintf('I,%d',i);
        end
        %integral features
        for i=1:4, fn{len+19+i}=sprintf('P,%d',i);
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
    
    if length(aux) > 10;
        %   warning off
        [hh,pp,lstat]=lillietest(aux);
        %    warning on
        sd(i)=lstat;
    else
        sd(i)=0;
    end
end
% fprintf('Lstat: ');
% fprintf('%1.2f ',sd);
% fprintf('\n');
warning('on','stats:lillietest:OutOfRangeP');
[max indA]=sort(-sd(1:len/2));
feature_namesA={fn{ind(1:inputs)}};
inspkA=cc(:,ind(1:inputs));
[max indB]=sort(-sd(len/2+1:len));
feature_namesB={fn{ind(1:inputs)}};
inspkB=cc(:,ind(1:inputs));


feature_namesC={fn{ind(1:inputs)}};
inspkC=cc(:,ind(1:inputs));
[max indD]=sort(-sd(len/2+1:len));
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

