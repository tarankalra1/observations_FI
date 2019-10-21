clear all ; close all ; clc; 
% TSK 
% Comparing three bedload formulations
% MPM with average waveform           
nt1=1; nt2= 2044; 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

% VSPEC 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght.mat',......
                                                           'bedldx_wh_vspec','R','Beta','Ur')
bedldx_empirical_vand_all=bedldx_wh_vspec; 
cumbedld_empirical_vand_all=cumtrapz(bedldx_empirical_vand_all); 

% WITHOUT SURFACE WAVE TERMS  that change wave period and without the reynolds stress term.  
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght_withoutwave10.mat',......
                                                           'bedldx_wh_vspec','R','Beta','Ur')
bedldx_empirical_vand_10=bedldx_wh_vspec; 
cumbedld_empirical_vand_10=cumtrapz(bedldx_empirical_vand_10);

% without surface wave terms that chagne wave period and reynolds stress term
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght_withoutwave11.mat',......
%                                                           'bedldx_wh_vspec','R','Beta','Ur')
%bedldx_empirical_vand_10=bedldx_wh_vspec; 
%cumbedld_empirical_vand_10=cumtrapz(bedldx_empirical_vand_10);

% without the change in wave period and without the reynolds stress term 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght_withoutwave00.mat',....
                                                                 'bedldx_wh_vspec','R','Beta','Ur')
bedldx_empirical_vand_00=bedldx_wh_vspec; 
cumbedld_empirical_vand_00=cumtrapz(bedldx_empirical_vand_00);

dt=3600 ; 
%  
bdld_plot=1; 
if(bdld_plot==1)
figure(1)
plot(dn(nt1:nt2),cumbedld_empirical_vand_all(nt1:nt2)*dt,'k'); 
% hold on 
% plot(dn(nt1:nt2),cumbedld_empirical_vand_10(nt1:nt2)*dt,'r'); 
hold on 
plot(dn(nt1:nt2),cumbedld_empirical_vand_00(nt1:nt2)*dt,'r'); 
%plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'r') 

legend('With surface wave streaming','Without surface wave streaming','Location','Northwest')
%legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('\int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig11_bedld_vspec_withoutwavestream.png')
end 
