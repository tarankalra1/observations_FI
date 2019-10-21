clear all ; close all ; clc; 
% TSK 
% Comparing three bedload formulations
% MPM with average waveform           
nt1=1; nt2= 2044; 
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

% VANDERA from workhorse using all the empircal waveform 

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur');
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat.mat',......
%    'bedldx_wh_vspec','R','Beta','Ur')

% % VSPEC 
% load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat_check.mat',....
%     'bedldx_wh_vspec','R','Beta','Ur')
% 
% cumbedld_empirical_vand=cumtrapz(bedldx_wh_vspec); 

% PSPEC - with phase lag + no surface stress terms
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat_pspec.mat',.....
%    'bedldx_wh_pspec','R','Beta','Ur')

% VSPEC 
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght.mat',......
                                                           'bedldx_wh_vspec','R','Beta','Ur')
%     
%  load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet_10cm.mat',.....
%     'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur')                                                   
%                                                        
%save('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat_pspec.mat',.....
%    'bedldx_wh_pspec','R','Beta','Ur')
bedldx_empirical_vand_all=bedldx_wh_vspec; 
cumbedld_empirical_vand_all=cumtrapz(bedldx_empirical_vand_all); 

% PSPEC WITHOUT SURFACE WAVE TERMS  
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat_pspec_minus_waveterms.mat',......
%    'bedldx_wh_pspec','R','Beta','Ur')
%bedldx_empirical_vand_minussur=bedldx_wh_pspec ; 
%cumbedld_empirical_vand_minussur=cumtrapz(bedldx_empirical_vand_minussur); 

 % PSPEC WITOUT PHASE LAG NO SURFACE TERMS 
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat_pspec_minus_waveterms_minusphase.mat',........
%    'bedldx_wh_pspec','R','Beta','Ur')
%bedldx_empirical_vand_minussur_minusphase=bedldx_wh_pspec ; 
%cumbedld_empirical_vand_minussur_minusphase=cumtrapz(bedldx_empirical_vand_minussur_minusphase); 
% 

% MPM bedload actual using full wave burst data
load('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_only_ss.mat'); 
cumbedld_mpm=cumtrapz(qb_measured);  
bedld_mpm=qb_measured; 

dt=3600 ; 
%  
bdld_plot=0; 
if(bdld_plot==1)
  figure(1)
subplot(2,1,1)
plot(dn(nt1:nt2),bedldx_wh_vspec(nt1:nt2),'k')
hold on
plot(dn(nt1:nt2),qb_measured(nt1:nt2),'r')
legend('Asymmetric bedload','Traditional bedload','Location','Best')
xlim([dn(nt1) dn(nt2)]); 
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('q_{b} (m^{2}/s)')
ylim([-2.0e-07 2.0e-05])
%figure(2) 

subplot(2,1,2)
plot(dn(nt1:nt2),cumbedld_empirical_vand_all(nt1:nt2)*dt,'k'); 
hold on 
plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'r') 
legend('Asymmetric bedload','Traditional bedload','Location','Northwest')
%legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
  datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
%ylabel('qb dt (m^{2})
ylabel('\int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
% \int_{0}^{2} x^2\sin(x) dx
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
 print('-dpng','-r300','fig5_bedld_vspec.png')

end 

%
%print -dpng 'fig5_cumsum_3600dt.png'

% % %     
% figure(2)
% plot(dn(nt1:nt2),bedldx_empirical_vand_minussur(nt1:nt2)-......
%                  bedldx_empirical_vand_minussur_minusphase(nt1:nt2))
%  xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) %            
