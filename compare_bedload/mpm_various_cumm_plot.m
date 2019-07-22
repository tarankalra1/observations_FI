clear all ; close all ; clc; 

% Load all the MPM files..

nt1=1; nt2=2044; 
 
load('../../matfiles/skewness_steve.mat','dn')
%IG+mean+incident
load('../../matfiles/mpm_only_ss_IG_addmean.mat','Su_skewness','qb_measured','bedload_measured'); 

bedload_measured_IG_nodtrend=bedload_measured;
qb_measured_IG_nodtrend=qb_measured; 
Su_measured_IG_nodtrend=Su_skewness; 

%IG +incident --after detrending 
load('../../matfiles/mpm_only_ss_IG.mat','Su_skewness','qb_measured','bedload_measured'); 
bedload_measured_IG=bedload_measured;
qb_measured_IG=qb_measured; 
Su_measured_IG=Su_skewness;

%incident --after detrending+deIG 
load('../../matfiles/mpm_only_ss.mat','Su_skewness','qb_measured','bedload_measured'); 
bedload_measured=bedload_measured;
Su_measured=Su_skewness; 

nt1=1; nt2=2044; 

dt=3600; 

figure(1)
plot(dn(nt1:nt2),bedload_measured(nt1:nt2)*dt,'k') 
hold on 
plot(dn(nt1:nt2),bedload_measured_IG(nt1:nt2)*dt,'b') 
hold on 
plot(dn(nt1:nt2),bedload_measured_IG_nodtrend(nt1:nt2)*dt,'r') 
hold on 
plot(dn(nt1:nt2),bedload_measured_IG_nodtrend(nt1:nt2)*dt*0.0,'r--') 
title('MPM bedload (m^{2})')
legend('only incident','incident+IG','incident+IG+mean','Location','Southwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
print('-dpng','-r300','../../pngfiles/cumulative_bedload_mpm_various.png')

figure(2)
plot(dn(nt1:nt2),qb_measured(nt1:nt2),'k') 
hold on 
plot(dn(nt1:nt2),qb_measured_IG(nt1:nt2),'b') 
hold on 
plot(dn(nt1:nt2),qb_measured_IG_nodtrend(nt1:nt2),'r') 
hold on 
plot(dn(nt1:nt2),qb_measured_IG_nodtrend(nt1:nt2),'r--') 
title('MPM bedload flux (m^{2}/s)')
legend('only incident','incident+IG','incident+IG+mean','Location','Southwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
%ylim([-0.002 0.008])
print('-dpng','-r300','../../pngfiles/bedload_flux_mpm_various.png')

figure(3)
plot(dn(nt1:nt2),Su_measured(nt1:nt2),'k') 
hold on 
plot(dn(nt1:nt2),Su_measured_IG(nt1:nt2),'b') 
hold on 
plot(dn(nt1:nt2),Su_measured_IG_nodtrend(nt1:nt2),'r') 
hold on
plot(dn(nt1:nt2),Su_measured_IG_nodtrend(nt1:nt2)*0.0,'r--') 
title('Su skewness)')
legend('only incident','incident+IG','incident+IG+mean','Location','Southwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
%ylim([-0.002 0.008])
print('-dpng','-r300','../../pngfiles/Su_skewness_various.png')
