clear all ; close all ; clc; 
% Comparing three bedload formulations
% MPM
% Direct waveform with averaged waveforms manually from ADV data
% all direct waveforms from ADV data
% STEVE's claculations of directly obtained ADV data 
nt1=1; nt2= 2044; 
% % 
% % 
% load('../../matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%             'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar',...
%             'ur_cube','ang_rot','Au_skewness')
%  
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 
     
 % STEVE's WAVEFORM
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_directwaveform.mat','bedldx','R','Beta','Ur');
bedld_directwaveform=bedldx; 
cumbedld_directwaveform=cumtrapz(bedldx);  

% VANDERA from workhorse using all the empircal waveform 

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur');
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_allwaveforms.mat','bedldx_allwaveform','R','Beta','Ur')
cumbedld_allwaveform=cumtrapz(bedldx_allwaveform);  

% MPM bedload 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_only_ss.mat'); 
cumbedld_mpm=cumtrapz(qb_measured);  
bedld_mpm=qb_measured; 

dt=3600 ; 

% 
figure(1) 
% plot(cumbedld_directwaveform(1:end),'k') 
 hold on
plot(cumbedld_allwaveform(1:end)*3600,'r'); 
%hold on 
%plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'b') 
legend('waveform averaged','all waveforms direct','Location','Northwest')
%legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
% print -dpng '../pngfiles/bedload/cumulative_bedload_vand_vspec_mpm_3600.png'
%print -dpng '../pngfiles/bedload/cumulative_bedload_vand_waveforms.png'

% % %   
% figure(2)
% plot(dn(nt1:nt2),bedld_directwaveform(nt1:nt2),'r')
% hold on 
% plot(dn(nt1:nt2),bedldx_allwaveform,'b')
% legend('waveform averaged','all waveforms direct','Location','Northeast')
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % 
% print -dpng '../pngfiles/bedload/bedload_vand_waveforms.png'
