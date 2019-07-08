clear z x ti
x=UBS(100).ur;
len=length(x);
fs=8;
idx=1
cnt=0
c=sign(x(idx))
while idx <  length(x)
cnt=cnt+1;
c=sign(x(idx));
while c==sign(x(idx)) & idx < length(x)
idx=idx+1;
end
z(cnt)=interp1([x(idx-1) x(idx)],[idx-1 idx],0);
end


ti=[1:len z(~isnan(z))];
[ti,I]=sort(ti);

x=[x' zeros(1,length(z(~isnan(z))))];
x=x(I);
t=ti/8;
plot(t,x)
nz=length(find(x==0))