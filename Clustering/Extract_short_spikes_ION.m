function [spikes,index]=Extract_short_spikes_ION(handles)

sr=handles.par.sr;
w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
int_factor=handles.par.int_factor;
% LOAD DATA
if exist(handles.filename,'file'),
   f=fopen(handles.filename,'r');
   x=fread(f,'int16');
   fclose(f);
   x=x(1:sr*60);

   x=x(:)';
   x=x*handles.par.transform_factor;%%in microvolts
end

switch handles.filterset
   case {-1,'1_5000'}, %no additional filters at all
      [b,a]=butter(2,[1 5000]/sr*2);
      x=filtfilt(b,a,x);
      fmin=15000;
   case {0,'none'}, %no additional filters at all filters at all
      fmin = 0; %add some basic filtering
   case {1,'unfilt1000',1000},
      fmin = 1000;
      x=b2_unfilt(x,sr,fmin,[]);
   case {2,'unfilt300',300},
      fmin = 300;
      x=b2_unfilt(x,sr,fmin,[]);
   case {3,'unfilt300filt1000'},
      fmin=300;
      x=b2_unfilt(x,sr,fmin,[]);
      [b,a]=butter(2,1000/sr*2,'high');
      x=filter(b,a,x);
      fmin=3001000;
   case {4,'unfilt300filt1000unfilt1000'},
      fmin=300;
      x=b2_unfilt(x,sr,fmin,[]);
      [b,a]=butter(2,1000/sr*2,'high');
      x=filter(b,a,x);
      fmin=1000;
      x=b2_unfilt(x,sr,fmin,[]);
      fmin=311000;
   case {5,'filt300'}, 
      fmin=300; %emulate causal filttering
      [b,a]=butter(2,fmin/sr*2,'high');
      x=filter(b,a,x);
      fmin=3000;
   case {6,'filt1000'}, 
      fmin=1000; %emulate causal filttering
      [b,a]=butter(2,fmin/sr*2,'high');
      x=filter(b,a,x);
      fmin=10000;
end
[b,a]=ellip(2,0.1,40,[100 3000]/sr*2);
x=filtfilt(b,a,x);

tens_fname=sprintf('tens_%s-%d.mat',handles.bname,handles.channel+fmin);
y=x;
x=y(1:min(round(sr*10),length(y)));
save(tens_fname,'x','sr');
x=y;
clear y

thr = handles.par.stdmin*median(abs(x))/0.6745;

% LOCATE SPIKE TIMES
switch handles.threshold
   case 'pos', ups = find(x(w_pre+2:end-w_post-2) > thr) +w_pre+1;%indices of values above threshold
   case 'neg', ups = find(x(w_pre+2:end-w_post-2) < -thr) +w_pre+1;%indices of values below threshold
end
if isempty(ups), %no spikes detected
   spikes=[]; index=[]; %#ok<NASGU>
   %    save(spikes_fname,'spikes','index');
   return;
end
%find beginnings and endings of intervals with values above threshold
begs=[0 find(diff(ups)>1)]+1;%index in index array ups
ends=[begs(2:end)-1 length(ups)];
begs=ups(begs);%indices in data array
ends=ups(ends);
%find two events which are closer than 0.2ms, remove smallest from two
%closest
bad=[];
for i=1:length(begs)-1,
   if begs(i+1)-ends(i)<=0.2e-3*sr, bad=[bad i]; end
end
nspk=length(begs);
index=NaN*ones(nspk,1);
np=w_post+w_pre;
spikes=NaN*ones(nspk,np+4);
m=NaN*ones(nspk,1);
for i=1:nspk,
   [m(i),ind] = max(abs(x(begs(i):ends(i))));
   index(i) = begs(i) + ind -1;
   spikes(i,:)= x (index(i)-w_pre+1-2:index(i)+w_post+2);
end
bad1=[];
for i=bad,
   if m(i)>m(i+1), bad1=[bad1 i+1];
   else bad1=[bad1 i]; end
end
ind=setdiff(1:nspk,bad1);
spikes=spikes(ind,:);
index=index(ind);
nspk=length(ind);
%INTERPOLATION
if handles.par.interpolation=='y',
   %Does interpolation
   ints=1/int_factor:1/int_factor:np+4;
   intspikes = spline(1:np+4,spikes,ints);
   spikes1=zeros(nspk,np);

   switch handles.threshold
      case 'pos'
         [maxi iaux]=max(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
      case 'neg'
         [maxi iaux]=min(intspikes(:,(w_pre+1)*int_factor:(w_pre+3)*int_factor),[],2);
   end
   iaux = iaux + (w_pre+1)*int_factor -1;
   for i=1:nspk
      spikes1(i,w_pre:-1:1) = intspikes(i,iaux(i):-int_factor:iaux(i)-(w_pre-1)*int_factor);
      spikes1(i,w_pre+1:np) = intspikes(i,iaux(i)+int_factor:int_factor:iaux(i)+w_post*int_factor);
   end
   spikes=spikes1;
   clear spikes1
else
   spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
   spikes(:,1:2)=[];
end
index=index/sr*1000;
% save(spikes_fname,'spikes','index');

   function d=b2_unfilt(d,sr,fmin,fmax)
      %filter signal backwards in time
      %butterworth 2nd order
      d=d(end:-1:1);
      if ~isempty(fmin)&&~isnan(fmin),
         [b,a]=butter(2,fmin/sr*2,'high');
         d=filter(b,a,d);
      end
      if ~isempty(fmax)&&~isnan(fmax),
         [b,a]=butter(2,fmax/sr*2,'low');
         d=filter(b,a,d);
      end
      d=d(end:-1:1);
   end

if 0
   %interpolation test for 
   w_pre=3;
   w_post=5;
   int_factor=2;
   np=w_post+w_pre;
   nspk=1;
   spikes=rand(nspk,np+4); 
   spikes=spikes+repmat([linspace(1,2.7,w_pre+2) linspace(1.5,1,w_post+2)],[nspk,1]);
   ints=1/int_factor:1/int_factor:np+4;
   intspikes = spline(1:np+4,spikes,ints);
   spikes1=zeros(nspk,np);

   
   
   [maxi iaux]=max(intspikes(:,(w_pre-2)*int_factor:(w_pre+2)*int_factor),[],2);
   iaux = iaux + (w_pre-2)*int_factor -1;
   for i=1:nspk
      spikes1(i,w_pre:-1:1) = intspikes(i,iaux(i):-int_factor:iaux(i)-(w_pre-1)*int_factor);
      spikes1(i,w_pre+1:np) = intspikes(i,iaux(i)+int_factor:int_factor:iaux(i)+w_post*int_factor);
   end
   i=1;
   clf;
   plot(1-2:np+2,spikes(i,:),'b');
   hold on; 
   plot(ints-2,intspikes(i,:),'g')

   plot(1:np,spikes1,'r');
end
end