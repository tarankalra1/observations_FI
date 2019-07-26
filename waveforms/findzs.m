function tz=findzs(s,t)
%function tz=findzs(s,t)
%
%m-fcn to find zero crossing times of a signal s
%
%INPUT:
% s = signal
% t = time
len=length(s);

idx=1
cnt=0
c=sign(s(idx))
while idx < len
cnt=cnt+1;
c=sign(s(idx));
while c==sign(s(idx)) & idx < len
idx=idx+1;
end
tz(cnt)=interp1([s(idx-1) s(idx)],[t(idx-1) t(idx)],0);
end


% ti=[1:len z(~isnan(z))];
% [ti,I]=sort(ti);
% 
% x=[x' zeros(1,length(z(~isnan(z))))];
% x=x(I);
% t=ti/8;
% plot(t,x)
% nz=length(find(x==0))