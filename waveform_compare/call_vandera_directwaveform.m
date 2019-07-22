clear all ; close all ; clc; 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform.mat',.....
     'Ur_emp','Hs','Td','h',.......
     'umax_emp','umin_emp','Tc_emp','Tt_emp',........
     'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp'); % 
% Read in Steve's waveform to get umax, umin, T_c, T_t....
%
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/9917adv_wfr.mat')

umax=[wfr.umax] ; 
umin=[wfr.umin] ;
T_c=[wfr.Tc]   ;
T_t=[wfr.Tt]   ;
T_cu=[wfr.Tcu] ;
T_tu=[wfr.Ttu] ; 
T=[wfr.T] ; 
R=[wfr.R] ;
uhat=[wfr.Uw] ;
%  
nt1=1; nt2=2044;

 for i=nt1:nt2
      if(~isnan(h(i)))
        Hs_emp(i)=Hs(i); 
        T_emp(i)=Td(i); 
        h_emp(i)=h(i); 
      end
 end
 figure(1)
 plot(dn(nt1:nt2),T_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T,'r')
title('Time period')
 legend('workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
print -dpng '../pngfiles/workhorse_waveform_ADV_Timeperiod.png'

figure(2)
 plot(dn(nt1:nt2),Tc_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_c,'r')
 title('Tcrest')
 legend('empirical from workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_T_crest.png'
 
 figure(3)
 plot(dn(nt1:nt2),uhat_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),uhat,'r')
 title('uhat')
 legend('empirical from workhorse','waveform from ADV')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_uhat.png'

 
 figure(4)
 plot(dn(nt1:nt2),umax_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),umax,'r')
 title('umax/ucrest')
 legend('empirical from workhorse','waveform from ADV')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_umax.png'
 
 
 figure(5)
 plot(dn(nt1:nt2),umin_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),umin,'r')
 title('umin/utrough')
 legend('empirical from workhorse','waveform from ADV','Location','Southeast')
 %xlim([dn(nt1) dn(nt2)]);
 %datetick('x',2) % 
 xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_umin.png'
 
 
 figure(6)
 plot(dn(nt1:nt2),RR_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),R,'r')
 title('R')
 legend('empirical from workhorse','waveform from ADV','Location','Southeast')
 %xlim([dn(nt1) dn(nt2)]);
 %datetick('x',2) % 
 xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_R.png'
 
 
figure(7)
 plot(dn(nt1:nt2),Tt_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_t,'r')
 title('T-trough')
 legend('empirical from workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_T_trough.png'
 
 
figure(8)
 plot(dn(nt1:nt2),Tcu_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_cu,'r')
 title('T-cu')
 legend('empirical from workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_T_cu.png'
 
figure(9)
 plot(dn(nt1:nt2),Ttu_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_tu,'r')
 title('T-tu')
 legend('empirical from workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/workhorse_waveform_ADV_T_tu.png'