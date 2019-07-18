function pltwaveform(wf,nb)
xlim=[0 16]

for i=1:length(wf)
plot(wf(i).t,wf(i).ub)
title(sprintf('Burst %d waveform %d: t0 @ %s',nb,i,datestr(wf(i).dn(1))))
set(gca,'xlim',xlim)
set(gca,'ylim',[-1.2 1.2])
xlabel('time, [s]')
ylabel('u_b, [m ^. s^{-1}]')
line(xlim,[0 0],'color',[0.5 0.5 0.5])
saveas(gcf,sprintf('9917adv_b%d_wf%d.png',nb,i),'png')
pause(0.3)

end

