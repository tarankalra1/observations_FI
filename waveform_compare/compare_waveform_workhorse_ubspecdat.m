clear all ; close all ; clc; 
% 
% WORKHORSE 
% Used ubspecdat. m to obtain empirical waveform with vspec and compare with Steve's manually
% obtained waveforms.
% Tarandeep S Kalra
% 

%    'WORKHORSE UBSPECDAT using vspec'
load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec.mat',.....
     'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
      'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');

 load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_wfr.mat')

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
      %if(~isnan((i)))
        Hs_emp(i)=Hs(i); 
        T_emp(i)=Tbr(i); 
        h_emp(i)=h(i); 
      %end
 end
figure(1)
plot(dn(nt1:nt2),T_emp(nt1:nt2),'k--')
hold on
plot(dn(nt1:nt2),T,'r')
title('Time period')
legend('empircal waveform workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_Timeperiod.png'

figure(2)
 plot(dn(nt1:nt2),Tc_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_c,'r')
 title('Tcrest')
 legend('empircal waveform workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_T_crest.png'
 
 figure(3)
 plot(dn(nt1:nt2),uhat_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),uhat ,'r')
 title('uhat')
 legend('empircal waveform workhorse','waveform from ADV')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_uhat.png'

 
 figure(4)
 plot(dn(nt1:nt2),umax_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),umax,'r')
 title('ucrest')
 legend('empircal waveform workhorse','waveform from ADV')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_umax.png'
 
 
 figure(5)
 plot(dn(nt1:nt2),umin_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),umin,'r')
 title('utrough')
 legend('empircal waveform workhorse','waveform from ADV','Location','Southeast')
 %xlim([dn(nt1) dn(nt2)]);
 %datetick('x',2) % 
 xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_umin.png'
 
%  
%  figure(6)
%  plot(dn(nt1:nt2),RR_emp(nt1:nt2),'k--')
%  hold on
%  plot(dn(nt1:nt2),R,'r')
%  title('R')
%  legend('empircal waveform workhorse','waveform from ADV')
%  %xlim([dn(nt1) dn(nt2)]);
%  %datetick('x',2) % 
%  xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % 
%  print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_R.png'
%  
 
figure(7)
 plot(dn(nt1:nt2),Tt_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_t,'r')
 title('T-trough')
 legend('empircal waveform workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_T_trough.png'
 
 
figure(8)
 plot(dn(nt1:nt2),Tcu_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_cu,'r')
 title('Tcu')
 legend('empircal waveform workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_T_cu.png'
 
figure(9)
 plot(dn(nt1:nt2),Ttu_emp(nt1:nt2),'k--')
 hold on
 plot(dn(nt1:nt2),T_tu,'r')
 title('T-tu')
 legend('empircal waveform workhorse','waveform from ADV')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 
 print -dpng '../pngfiles/waveform_compare/waveform_whvspectra_T_tu.png'
