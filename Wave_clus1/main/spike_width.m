function [w1, w2, amp, halfamp]=spike_width(m)
%define width parameters as in Bartho et al 2004
%w1 width from min to max
%w2 width from max to min

warning off MATLAB:m_warning_end_without_block

n=length(m);
factor=10;
time=1:n;
time1=1:1/factor:n;
m=interp1(time,m,time1,'spline');


[mn,indmn]=min(m);
[mx,indmx]=max(m(indmn+1:end));
indmx=indmx+indmn;

w1=(indmx-indmn)/factor;
amp=mx-mn;

[mx,indmx]=max(m(1:indmn));
w2=(-indmx+indmn)/factor;

halfamp=0;

% bl=mean(m(1:floor(indmn/2)));
% halfamp=(mn-bl)/2;
% 
% indhalf=find(m(1:indmx)<halfamp);
% w2=length(indhalf)/factor;

end

% old version 
% %define width parameters as in Bartho et al 2004
% %w1 width from min to max
% %w2 half amplitude duration
% 
% warning off MATLAB:m_warning_end_without_block
% 
% n=length(m);
% factor=10;
% time=1:n;
% time1=1:1/factor:n;
% m=interp1([time],[m],time1,'spline');
% 
% 
% [mn,indmn]=min(m);
% [mx,indmx]=max(m(indmn+1:end));
% indmx=indmx+indmn;
% 
% w1=(indmx-indmn)/factor;
% amp=mx-mn;
% 
% bl=mean(m(1:floor(indmn/2)));
% halfamp=(mn-bl)/2;
% 
% indhalf=find(m(1:indmx)<halfamp);
% w2=length(indhalf)/factor;

