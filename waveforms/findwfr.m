function wfr=findwfr(wf)

% for i=1:length(wf)
% ubh(i)=wf(i).umax-wf(i).umin;
% T(i)=wf(i).T;
% end
% Uw=mean(ubh)/2;
% T=mean(T);

for i=1:length(wf)
Uw(i)=(wf(i).umax-wf(i).umin)/2;
T(i)=wf(i).T;
Tr(i)=wf(i).T/T(i);    
Tcr(i)=wf(i).Tc/T(i);
Ttr(i)=wf(i).Tt/T(i);
Acr(i)=wf(i).Ac/(Uw(i)*T(i));
Atr(i)=wf(i).At/(Uw(i)*T(i));
umaxr(i)=wf(i).umax/Uw(i);
uminr(i)=wf(i).umin/Uw(i);
Tcur(i)=wf(i).Tcu/T(i);
Ttur(i)=wf(i).Ttu/T(i);
end

wfr.dn=wf(1).dn(1);
wfr.Uw=mean(Uw);
wfr.T=mean(T);
wfr.umax=mean(umaxr)*wfr.Uw;
wfr.umin=mean(uminr)*wfr.Uw;
wfr.Tc=mean(Tcr)*wfr.T;
wfr.Tt=mean(Ttr)*wfr.T;
wfr.Tcu=mean(Tcur)*wfr.T;
wfr.Ttu=mean(Ttur)*wfr.T;
wfr.Ac=mean(Acr)*(wfr.Uw*wfr.T);
wfr.At=mean(Atr)*(wfr.Uw*wfr.T);
wfr.R=wfr.umax/(wfr.umax-wfr.umin);
wfr.alpha=2*wfr.Tcu/wfr.T;