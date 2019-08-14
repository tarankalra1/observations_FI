clear all ; close all ; clc; 
% 
% WORKHORSE 
% Used ubspecdat. m to obtain empirical waveform with vspec and compare with Steve's manually
% obtained waveforms.
% Tarandeep S Kalra
% 
% Fidn the reason for difference 

%    'WORKHORSE UBSPECDAT using vspec' 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_pspec_inithght.mat',.....
      'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
       'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp'); 
   
% load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') ; 

 % ADV calculated major angle 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','dn','ang_rot');

load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_wfr.mat'); 

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
beginTime='07-02-2014'; finishTime='05-05-2014';
date_frmt='dd-mm-yyyy';

figure(1)
%subplot(3,1,1)
plot(dn(nt1:nt2),T(nt1:nt2)-T_emp(nt1:nt2),'k--')
%title('del T ')
% legend('empircal waveform workhorse','waveform from ADV')
%xlim([dn(nt1) dn(nt2)]);datetick('x',2) % 
hold on 

%subplot(3,1,2)
plot(dn(nt1:nt2),uhat(nt1:nt2),'r')
datetick('x',2)
set(gca,'XLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)])
%xlim([dn(nt1) dn(nt2)]);datetick('x',2) % 
  legend('del T', 'Hs')
% 
% figure(2)
% %subplot(3,1,1)
% plot(dn(nt1:nt2),uhat(nt1:nt2)-uhat_emp(nt1:nt2),'k--')
% %title('del T ')
% % legend('empircal waveform workhorse','waveform from ADV')
% %xlim([dn(nt1) dn(nt2)]);datetick('x',2) % 
% hold on 
% 
% %subplot(3,1,2)
% plot(dn(nt1:nt2),Hs(nt1:nt2),'r')
% datetick('x',2)
% set(gca,'XLim',[datenum(beginTime,date_frmt),datenum(finishTime,date_frmt)])
% %xlim([dn(nt1) dn(nt2)]);datetick('x',2) % 
% legend('del T', 'Hs')

%subplot(3,1,3)
%plot(dn(nt1:nt2),ang_rot(nt1:nt2),'k')
%title('Angle along major wave dir')
%xlim([dn(nt1) dn(nt2)]);datetick('x',2) % 

% print -dpng '../pngfiles/waveform_compare/waveform_whpspectra_T_crest.png'


