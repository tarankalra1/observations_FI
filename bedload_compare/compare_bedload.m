clear all ; close all ; clc; 
% Comparing three bedload formulations
% MPM, vandera with ubspecfun (JONSWAP), vandera with averaged waveformfrom
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

load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur');
%load('vandera_bedld_workhorse_ubspecdat.mat','bedldx_wh_vspec','R','Beta','Ur')

bedld_empirical_vand=bedldx; 
cumbedld_empirical_vand=cumtrapz(bedld_empirical_vand);  

% MPM bedload 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_only_ss.mat'); 
cumbedld_mpm=cumtrapz(qb_measured);  
bedld_mpm=qb_measured; 

dt=3600 ; 

% 
% bedload_mpm=cumtrapz(
% % % % % integrate in time 
%  load('../../matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur')
% subplot(2,1,2)
%  bedload_vand=cumtrapz(bedldx);
figure(1) 
plot(dn(nt1:nt2),cumbedld_directwaveform(nt1:nt2)*dt,'k') 
hold on
plot(dn(nt1:nt2),cumbedld_empirical_vand(nt1:nt2)*dt,'r'); 
hold on 
plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'b') 
legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/bedload/cumulative_bedload_vand_mpm_3600.png'

% %   
figure(2)
% %subplot(2,1,1)
plot(dn(nt1:nt2),bedld_directwaveform(nt1:nt2),'k')
hold on
plot(dn(nt1:nt2),bedld_empirical_vand(nt1:nt2),'r')
hold on 
plot(dn(nt1:nt2),bedld_mpm,'b')
legend('direct waveform','emp. waveform vand','mpm')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/bedload/bedload_vand_mpm.png'

figure(3) 
plot(dn(nt1:nt2),cumbedld_directwaveform(nt1:nt2)*dt,'r'); 
hold on 
plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'b') 
legend('direct waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/bedload/bedload_vand_direct_mpm.png'


% % title('Cumulative measured bedload with MPM')
%  xlim([dn(nt1) dn(nt2)]);
%  datetick('x',2) % 
% % 
% % % % integrate in time 
% load('../../matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur')
% hold on
% bedload_vand=cumtrapz(bedldx) ;
%  plot(dn(nt1:nt2),bedload_vand(nt1:nt2)*dt) 
% % title('Cumulative calculated with vandera')
%  xlim([dn(nt1) dn(nt2)]);
%  datetick('x',2) % 
%  legend('measured-MPM','vandera','location','NorthWest')
%  print -dpng '../pngfiles/cumulative_bedload_vand_mpm_withsign_3600.png'
% 
%  
