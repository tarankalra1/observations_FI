clear
load('mat\puv_proc_FI_iwaves.mat')
for ii=1:length(UBS)
b(ii).wf=findurwaveform(UBS(ii).ur',fs,dn(ii))
end

save mat\9917adv_waveforms b