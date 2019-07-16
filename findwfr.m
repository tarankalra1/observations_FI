function wfr=findwfr(wf)

for i=1:length(wf)
ubh(i)=wf(i).umax-wf(i).umin;
end
ubhr=mean(ubh)

for i=1:length(wf)
Tcr(i)=wf(i).Tc/ubhr;
Ttr(i)=wf(i).Tt/ubhr;
Acr(i)=wf(i).Ac/ubhr;
Atr(i)=wf(i).At/ubhr;
umaxr(i)=wf(i).umax/ubhr;
uminr(i)=wf(i).umin/ubhr;
tmax(i)=wf(i).tc_max/wf(i).Tc;
tmin(i)=wf(i).tt_max/wf(i).Tt;
end

wfr.dn=wf(1).dn(1);
wfr.ubh=ubhr;
wfr.umax=mean(umaxr)*ubhr;
wfr.umin=mean(uminr)*ubhr;
wfr.Tc=mean(Tcr)*ubhr;
wfr.Tt=mean(Ttr)*ubhr;
wfr.Ac=mean(Acr)*ubhr;
wfr.At=mean(Atr)*ubhr;
wfr.tmax=mean(tmax)*wfr.Tc;
wfr.tmin=mean(tmin)*wfr.Tt;
wfr.R=wfr.umax/(wfr.umax-wfr.umin);