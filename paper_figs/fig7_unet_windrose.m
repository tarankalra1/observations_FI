clear all ; close all ; clc; 
% Tarandeep S Kalra 
% Code to plot the crest/trough velocity/ wind rose of currents 

% Get the crest and trough velocity
load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec_inithght.mat',.....
      'umax_emp','umin_emp');  
   
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

% get the current velocity magnitude
load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/unet_9917adv_phi.mat',.......
                                               'phid','wdir','cdir','unet')

umag_curr=unet  ; % unet;  %.2814;% abs(0.0);
phi_curwave=phid;% 79.92*deg2rad ;% 0.0*deg2rad;

nt1=1 ; nt2=2044; 

figure(1)
%plot(dn(nt1:nt2),umax_emp(nt1:nt2),'k');
%hold on
%plot(dn(nt1:nt2),umin_emp(nt1:nt2),'r');
%hold on
subplot(2,1,1)
plot(dn(nt1:nt2),unet(nt1:nt2),'k');
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
%legend('Crest','Trough','Current')
%xlabel('Time (year/month/day)')
ylabel('Current velocity (m/s)')
%print -dpng '../pngfiles/repres/uhat_adv_wh_vspecdat.png'
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
%print('-dpng','-r300','fig1_ubr.png')
%figure(2)
subplot(2,1,2)
plot(dn(nt1:nt2), phid(nt1:nt2)) % ,'k','filled')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
%legend('Crest','Trough','Current')
xlabel('Time (year/month/day)')
ylabel('{\phi} (Angle between waves and currents)')
%print -dpng '../pngfiles/repres/uhat_adv_wh_vspecdat.png'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
%print('-dpng','-r300','fig7b_unet_phidir.png')

reg_windrose=0; 
if(reg_windrose==1)
 figure(2)
 WindRose(phi_curwave,umag_curr,'vwinds',[0 0.05 0.1 0.15],'TitleString',{''},.......
    'LabLegend','Current velocity (m/s)','LegendVariable','u_{d}' ,.....
    'labels',{'(90°)','(270°)','(0°)','(180°)'})
   print('-dpng','fig7_phidwave_2.png','-painters')
end 

scatter_plot=0;
if(scatter_plot==1)
figure(3)
WindRose(cdir,umag_curr,'vwinds',[0 0.05 0.1 0.15],'TitleString',{''},.......
    'LabLegend','Current velocity (m/s)','LegendVariable','u_{d}',.....
    'labels',{'North(0°)','South (180°)','East (90°)','West (270°)'})
print('-dpng','fig7_cdir.png','-painters')

figure(4)
ScatterWindRose(phi_curwave,umag_curr)
print('-dpng','-r300','fig7_phi_curwave_scatter.png')
 
figure(5)
ScatterWindRose(cdir,umag_curr)
print('-dpng','-r300','fig7_cdir_scatter.png')
end 
