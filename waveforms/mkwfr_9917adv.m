clear
load('mat\9917adv_waveforms')
for ii=1:length(b)
    wf=b(ii).wf;
    wfr(ii)=findwfr(wf);
end

save mat\9917adv_wfr wfr