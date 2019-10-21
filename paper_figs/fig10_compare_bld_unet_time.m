clear all ; close all ; clc; 
% TSK 
% Comparing the current + wave driven bedload formula 
%
% Get original time period for each crest and trough cycle 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec_inithght.mat',.....
      'Tc_emp','Tt_emp','RR_emp');

nt1=1; nt2= 2044; 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet.mat',.....
    'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc','Tt')
bedld_unet_1=bedldx_wh_unet; 
cumbedld_unet_1=cumtrapz(bedld_unet_1); 
Tc_final1=Tc; 
Tt_final1=Tt;
% Get the bedlaod and final crest and trough wave period
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet_timezero.mat',.....
%    'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc','Tt')
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet_timezero.mat',.....
     'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc_ini','Tc','Tt_ini','Tt')

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet_timezero.mat',.....
%     'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc_ini','Tc','Tt_ini','Tt') 
 
bedld_unet_2=bedldx_wh_unet; 
cumbedld_unet_2=cumtrapz(bedld_unet_2); 
Tc_final2=Tc; 
Tt_final2=Tt; 

dt=3600 ; 
%  
bdld_plot=1; 
if(bdld_plot==1)
figure(1)
plot(dn(nt1:nt2),cumbedld_unet_1(nt1:nt2)*dt,'k')
hold on
plot(dn(nt1:nt2),cumbedld_unet_2(nt1:nt2)*dt,'r')
legend('With half cycle modification','Without half cycle modification','Location','Best')
xlim([dn(nt1) dn(nt2)]); 
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('\int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig10_bedld_cum_unet_time.png') ; 
%ylim([-2.0e-07 2.0e-05])
%figure(2) 
end

time=1; 
if(time==1)
%subplot(2,1,2)
figure(2)
subplot(2,1,1)
plot(dn(nt1:nt2),Tc_final1(nt1:nt2),'k','LineWidth',2); 
hold on 
plot(dn(nt1:nt2),Tc_final2(nt1:nt2),'r') 
legend('With half cycle modification','Without half cycle modification','Location','Best')
ylabel('Crest cycle (s)')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')

subplot(2,1,2)
plot(dn(nt1:nt2),Tt_final1(nt1:nt2),'k','LineWidth',2); 
hold on 
plot(dn(nt1:nt2),Tt_final2(nt1:nt2),'r') 
legend('With half cycle modification','Without half cycle modification','Location','Best')
ylabel('Trough cycle (s)')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
 print('-dpng','-r300','fig10_crest_trough.png')
% 
% end 
% 
% figure(3)
% plot(Tc_emp+Tt_emp,'ko')
% hold on
% plot(Tc_final+Tt_final,'r')
end 