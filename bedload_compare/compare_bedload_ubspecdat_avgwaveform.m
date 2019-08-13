clear all ; close all ; clc; 
% Comparing three bedload formulations
% MPM with average waveform                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            , vandera with ubspecfun (UBSPECDAT),  
nt1=1; nt2= 2044; 
% %  ness')
%  
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 
 

% VANDERA from workhorse using all the empircal waveform 

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur');
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_ubspecdat.mat',......
    'bedldx_wh_vspec','R','Beta','Ur')

bedld_empirical_vand=bedldx_wh_vspec; 
cumbedld_empirical_vand=cumtrapz(bedld_empirical_vand);  

% MPM bedload 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_avgwaveform.mat',.......
                        'qb_measured','bedload_measured');
cumbedld_mpm=cumtrapz(qb_measured);  
bedld_mpm=qb_measured; 

dt=3600 ; 

%  
figure(1) 
% plot(dn(nt1:nt2),cumbedld_directwaveform(nt1:nt2)*dt,'k') 
% hold on
plot(dn(nt1:nt2),cumbedld_empirical_vand(nt1:nt2)*dt,'r'); 
hold on 
plot(dn(nt1:nt2),cumbedld_mpm(nt1:nt2)*dt,'b') 
legend('emp. waveform vand','mpm','Location','Northwest')
%legend('direct waveform vand','emp. waveform vand','mpm','Location','Northwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/bedload/cumulative_vand_vspec_mpm_avgwaveform_3600.png'

% %     