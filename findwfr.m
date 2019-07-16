function wfr=findwfr(wf)

for i=1:length(wf)
ubh(i)=wf(i).umax-wf(i).umin;
end
Uw=mean(ubh)/2;

for i=1:length(wf)
Tr(i)=wf(i).T/Uw;    
Tcr(i)=wf(i).Tc/Uw;
Ttr(i)=wf(i).Tt/Uw;
Acr(i)=wf(i).Ac/Uw;
Atr(i)=wf(i).At/Uw;
umaxr(i)=wf(i).umax/Uw;
uminr(i)=wf(i).umin/Uw;
Tcur(i)=wf(i).Tcu/wf(i).Tc;
Ttur(i)=wf(i).Ttu/wf(i).Tt;
end

wfr.dn=wf(1).dn(1);
wfr.Uw=Uw;
wfr.umax=mean(umaxr)*Uw;
wfr.umin=mean(uminr)*Uw;
wfr.T=mean(Tr)*Uw;
wfr.Tc=mean(Tcr)*Uw;
wfr.Tt=mean(Ttr)*Uw;
wfr.Tcu=mean(Tcur)*wfr.Tc;
wfr.Ttu=mean(Ttur)*wfr.Tt;
wfr.Ac=mean(Acr)*Uw;
wfr.At=mean(Atr)*Uw;
wfr.R=wfr.umax/(wfr.umax-wfr.umin);
wfr.alpha=2*wfr.Tcu/wfr.T;