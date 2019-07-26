function wf=stokes2ubwf(Hs,T,h)
%function wf=stokes2ubwf(Hs,T,h)
%generate stokes nearbed velocity wave form

t=0:T/128:2*T %generate 2 wave cycles
om=2*pi/T;
kh=qkhfs(om,h)
k=kh/h;
a=Hs/2;

%from http://pages.bangor.ac.uk/~oss062/RCEM%20Short%20Course/Davies/Excersices/Ex6_SedTr_Waves.pdf
U1=a*om*1/sinh(kh);
W=3*k*a/(4*sinh(kh)^3); %asymmetry parameter
ft=cos(om.*t)+W*cos(2*om.*t)
ub=U1*ft;

%find waveform wf
try
    %find zero crossing times
    tz=findzs(ub,t);


    %loop to find crest/trough info
    for ii=1:length(tz)-2
        %find half wave signal
        idx=find(t>=tz(ii) & t<=tz(ii+1));
        s(ii)={[0 ub(idx) 0]};
        ts(ii)={[tz(ii) t(idx) tz(ii+1)]};
        A(ii)=trapz(ts{ii},s{ii});
        S_max(ii)=max(s{ii});
        S_min(ii)=min(s{ii});
        T(ii)=tz(ii+1)-tz(ii);
    end

    cnt=0
    xlim=[0 15]
    if sign(A(1))>=0
        for i=1:2:length(ts)-1
        cnt=cnt+1
        %wf(cnt).dn=dnst+[ts{i} ts{i+1}(2:end)]/(3600*24);
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)]
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).ubmax=S_max(i);
        wf(cnt).ubmin=S_min(i+1);
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).ubmax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).ubmin)-wf(cnt).Tc;
    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end

    else
        for i=2:2:length(ts)-1
        cnt=cnt+1
        %wf(cnt).dn=dnst+[ts{i} ts{i+1}(2:end)]/(3600*24);
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)]
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).ubmax=S_max(i);
        wf(cnt).ubmin=S_min(i+1);
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).ubmax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).ubmin)-wf(cnt).Tc;

    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end
    end
catch
        %wf.dn=dnst;
        wf.t=NaN;
        wf.ub=NaN;
        wf.T=NaN;
        wf.Tc=NaN;
        wf.Tt=NaN;
        wf.Nmax=NaN;
        wf.Nmin=NaN;
        wf.Ac=NaN;
        wf.At=NaN;
        wf.Tcu=NaN;
        wf.Ttu=NaN;
end

return

%%
function kh = qkhfs( w, h )
% QKHFS - Quick iterative calculation of kh in dispersion relationship
% kh = qkhf( w, h )
%
% Input:
%  w Angular wave frequency = 2*pi/T where T = wave period [1/s]
%  h Water depth [m]
% Returns:
%  kh = wavenumber * depth [ ]
% 
% Either w or h can be a vector, but not both.
% Hard-wired for MKS units.
% Orbital velocities from kh are accurate to 3e-12 !
%
% RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
% HR Wallingford Report TR 155, February 2006
% Eqns. 12a - 14

% csherwood@usgs.gov
% Sept 10, 2006

g = 9.81;
x = w.^2*h./g;
y = sqrt(x) .* (x<1) + x.* (x>=1);
%this appalling bit of code is about 25% faster than a for loop
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
kh=y;
return;

%%
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