clear all ; close all ; clc; 
% written by: Tarandeep S Kalra 
% Compare the ADV based derived waveform and ADV based empirical waveform. 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/adv_emp_waveform.mat',.....
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
RR=[wfr.R] ;
uhat=[wfr.Uw] ;
%  
nt1=1; nt2=2044;

 for i=nt1:nt2
      %if(~isnan((i)))
        Hs_emp(i)=Hs(i); 
        T_emp(i)=Td(i); 
        h_emp(i)=h(i); 
      %end
 end

figure(1)
scatter(T,T_emp,'.')
%R=corrcoef(T,T_emp)
%RR=R(1,2)
%title(['Period, correlation coeff is ',num2str(RR)])
title('period')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_period.png'

figure(2)
scatter(umax,umax_emp,'.')
%R=corrcoef(umax,umax_emp)
%RR=R(1,2)
%title(['Ucrest, correlation coeff is ',num2str(RR)])
title('ucrest')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_ucrest.png'
 
 figure(3)
scatter(umin,umin_emp,'.')
%R=corrcoef(umin,umin_emp)
%RR=R(1,2)
%title(['Utrough, correlation coeff is ',num2str(RR)])
title('Utrough')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_utrough.png'

figure(4)
scatter(T_c,Tc_emp,'.')
%R=corrcoef(T_c,Tc_emp)
%RR=R(1,2)
%title(['Tcrest, correlation coeff is ',num2str(RR)])
title('Tcrest')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_tcrest.png'
 
figure(5)
scatter(T_t,Tt_emp,'.')
%R=corrcoef(T_t,Tt_emp)
%RR=R(1,2)
%title(['Ttrough, correlation coeff is ',num2str(RR)])
title('Ttrough')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_ttrough.png'
 
figure(6)
scatter(RR,RR_emp,'.')
title('velocity skewness')
%R1=corrcoef(R,RR_emp)
%RR=R1(1,2)
%title(['Ttrough, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_Rskew.png'

figure(7)
scatter(T_cu,Tcu_emp,'.')
%R1=corrcoef(T_cu,Tcu_emp)
%RR=R1(1,2)
%title(['Tcu, correlation coeff is ',num2str(RR)])
title('Tcu')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_tcu.png'

figure(8)
scatter(T_tu,Ttu_emp,'.')
%R1=corrcoef(T_tu,Ttu_emp)
%RR=R1(1,2)
%title(['Ttu, correlation coeff is ',num2str(RR)])
title('Ttu')
xlabel('Direct-ADV')
ylabel('Empirical-ADV')
print -dpng '../pngfiles/waveform_compare/waveformscatter_ADV_ttu.png'
