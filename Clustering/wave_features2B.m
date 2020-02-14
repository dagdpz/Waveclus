function [inspk,feature_names,spikesadd,spikecomp] = wave_features2B(spikes,handles)
%Calculates the spike features

scales = handles.par.scales;
feature = handles.par.features;
inputs = handles.par.inputs;
int_factor = handles.par.int_factor;
w_pre = handles.par.w_pre;
w_post = handles.par.w_post;
sr = handles.par.sr;

% sr = 30000;
% handles.par.wavelet='haar'; 
% scales = 4;
% feature = 'wavpca';   
% inputs = 20;
% int_factor = 2;
% w_pre = 20;
% w_post = 44;
shift = 2;
shift = shift*int_factor;

nspk=size(spikes,1);
% len = size(spikes,2);
len = size(spikes,2)/int_factor;

offstart = w_pre * (int_factor -1)+shift;
offend  = w_post * (int_factor -1)+shift;
offend2 = offend + w_pre * int_factor;
downfac = 4;

%set(handles.file_name,'string','Calculating spike features ...');

% CALCULATES FEATURES
switch feature
    case 'wav'
        cc=zeros(nspk,len);
        for i=1:nspk                                % Wavelet decomposition
            [c,l]=wavedec(spikes(i,offstart+1:offend2),scales,handles.par.wavelet);
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
        %real spikes
        spikes1=zeros(nspk,len/(downfac/2));
        spikes1(:,(w_pre*int_factor -offstart)/(downfac/2):-1:1) = spikes(:,w_pre*int_factor:-(downfac/2):offstart+1);
        spikes1(:,((w_pre*int_factor)/(downfac/2)+1:offend2/(downfac/2))-offstart/(downfac/2)) = spikes(:,w_pre*int_factor+(downfac/2):(downfac/2):offend2);
        
        %wavelet
        cc=zeros(nspk,len/(downfac/2));
        for i=1:nspk                                % Wavelet decomposition
            [c,l]=wavedec(spikes1(i,:),scales,handles.par.wavelet);
            cc(i,1:len/(downfac/2))=c(1:len/(downfac/2));             %def.cc(i,1:ls)=c(1:ls);
        end
        spikecomp = cc(:,1:len/(downfac/2));
        clear spikes1
        
        %real spikes
        spikes1=zeros(nspk,len/downfac);
        spikes1(:,(w_pre*int_factor -offstart)/downfac:-1:1) = spikes(:,w_pre*int_factor:-downfac:offstart+1);
        spikes1(:,((w_pre*int_factor)/downfac+1:offend2/downfac)-offstart/downfac) = spikes(:,w_pre*int_factor+downfac:downfac:offend2);
        
        %         cc(1:nspk,len+1:2*len)=spikes(1:nspk,1:len);
        cc = cat(2,cc,spikes1);
        clear spikes1
        
        %pca
        [C,S] = princomp(spikes(:,offstart+1:offend2));
        cc = cat(2,cc,S(:,1:3));
        spikesadd = S(:,1:3);
        clear C S
        
%         %peaktopeak
%         6 18 27
%         [b,a] = butter(2,3000/(sr*int_factor/2),'low');
%         first = filtfilt(b,a,spikes(:,w_pre*int_factor:end)');
%         firstB = diff(first(5:end,:))<=0;
%         firstB = logical(cat(1,firstB,ones(1,size(first,2))));
%         peakindfirst = zeros(size(first,2),1);
%         firstamp = zeros(size(first,2),1);
%         for j = 1 : size(first,2)
%             dummy = 5:size(first,1);
%             peakindfirst(j) = min(dummy(firstB(:,j)));
%             clear dummy
%             firstamp(j) = first(peakindfirst(j),j)-spikes(j,20)';
%         end
%         
%         clear first firstB b a 
%         cc = cat(2,cc,peakindfirst,firstamp);
%         clear peakindfirst firstamp 
        %assign names of features
        %       fn=cell(1,len+6);
        
        fn=cell(1,len/(downfac/2)+len/downfac+3);
        
        %wavelet features
        j=1;
        k=scales+1;
        ll=cumsum(l);
        for i=ll(1:end-1),
            j1=1;
            while j<=i,
                if i==ll(1), fn{j}=sprintf('A%d,%d',k-1,j1);
                else fn{j}=sprintf('D%d,%d',k,j1);
                end
                j=j+1; j1=j1+1;
            end
            k=k-1;
        end
        %real spikes
        for i=1:len/downfac, fn{len/(downfac/2)+i}=sprintf('S,%d',i);
        end
        
        %pca features
        for i=1:3, fn{len/(downfac/2)+len/downfac+i}=sprintf('PCA,%d',i);
        end
%         %peaktopeak
%         for i=1:2, fn{2*len+3+i}=sprintf('PtP,%d',i);
%         end
        
end
warning('off','stats:lillietest:OutOfRangeP');
warning('off','stats:lillietest:OutOfRangePHigh');
warning('off','stats:lillietest:OutOfRangePLow');
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
warning('on','stats:lillietest:OutOfRangePHigh');
warning('on','stats:lillietest:OutOfRangePLow');

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

