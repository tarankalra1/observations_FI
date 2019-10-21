clear all ; close all ; clc; 
% TSK 
% Comparing the current + wave driven bedload formula 
% 
nt1=1; nt2= 2044; 
   load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

% VSPEC 
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght.mat',......
                                                           'bedldx_wh_vspec','R','Beta','Ur')
 % only wave driven part
bedldx_empirical_vand_all=bedldx_wh_vspec; 
cumbedld_empirical_vand_all=cumtrapz(bedldx_empirical_vand_all); 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet.mat',.....
    'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur')
bedld_unet_time=bedldx_wh_unet; 
cumbedld_bedld_unet_time=cumtrapz(bedld_unet_time); 


dt=3600 ; 
%  
bdld_plot=1; 
if(bdld_plot==1)
figure(1)
% subplot(2,1,1)
plot(dn(nt1:nt2),bedldx_empirical_vand_all(nt1:nt2),'k','LineWidth',2)
hold on
plot(dn(nt1:nt2),bedld_unet_time(nt1:nt2),'r--')
legend('Wave driven','Wave and current driven','Location','Best')
xlim([dn(nt1+1200) dn(nt2)]); 
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('q_{b} (m^{2}/s)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
%print('-dpng','-r300','fig6_bedld_cum_unet_1.png') ; 
%ylim([-2.0e-07 2.0e-05])
%figure(2) 

%subplot(2,1,2)
figure(2)
plot(dn(nt1:nt2),cumbedld_empirical_vand_all(nt1:nt2)*dt,'k'); 
hold on 
plot(dn(nt1:nt2),cumbedld_bedld_unet_time(nt1:nt2)*dt,'r') 
legend('Wave driven','Wave and current driven','Location','Northwest')
%legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1+1200) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
% %ylabel('qb dt (m^{2})
ylabel('\int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
% % \int_{0}^{2} x^2\sin(x) dx
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
% print('-dpng','-r300','fig6_bedld_cum_unet_2.png')

end 
