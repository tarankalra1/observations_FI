clear all ; close all ; clc; 
%
% WORKHORSE
%
% Used ubspecdat. m to obtain empirical waveform with vspec and compare with Steve's manually
% obtained waveforms.
% VSPEC
load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec.mat',.....
     'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
      'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform.mat',.....
%     'Ur_emp','Hs','Tbav','h',.......
%     'umax_emp','umin_emp','Tc_emp','Tt_emp',........
%     'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');

load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_wfr.mat')
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

 isave=1 ; 
if(isave==1)
% figure(1)
% scatter(T,T_emp,'.')
% %R=corrcoef(T,T_emp);
% %RR=R(1,2)
% title('Period')
% xlabel('Direct-ADV')
% ylabel('Empirical-WH')
% print -dpng '../pngfiles/waveform_compare/waveformscatter_wh_period.png'

figure(2)
scatter(umax,umax_emp,'.')
R=corrcoef(umax,umax_emp,'row','complete');
RR=R(1,2)
title(['ucrest, correlation coeff is ',num2str(RR)])
%title('Ucrest')
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_ucrest.png'
 
 figure(3)
scatter(umin,umin_emp,'.')
R=corrcoef(umin,umin_emp,'row','complete');
RR=R(1,2)
title(['utrough, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_utrough.png'

figure(4)
scatter(uhat,uhat_emp,'.')
R=corrcoef(uhat,uhat_emp,'row','complete');
RR=R(1,2)
title(['uhat, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_uhat.png'

figure(5)
scatter(T_c,Tc_emp,'.')
R=corrcoef(T_c,Tc_emp,'row','complete');
RR=R(1,2)
title(['Tcrest, correlation coeff is ',num2str(RR)])
%title('Tcrest')
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_tcrest.png'
 
figure(6)
scatter(T_t,Tt_emp,'.')
R=corrcoef(T_t,Tt_emp,'row','complete');
RR=R(1,2)
title(['Ttrough, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_ttrough.png'
 
figure(7)
scatter(T_cu,Tcu_emp,'.')
R=corrcoef(T_cu,Tcu_emp,'row','complete');
RR=R(1,2)
title(['Tcu, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_tcu.png'

figure(8)
scatter(T_tu,Ttu_emp,'.')
%scatter(T_tu,Ttu_emp,'.')
R=corrcoef(T_tu,Ttu_emp,'row','complete');
RR=R(1,2)
title(['Ttu, correlation coeff is ',num2str(RR)])
xlabel('Direct-ADV')
ylabel('Empirical-WH')
print -dpng '../pngfiles/waveform_compare/waveformscatter_whvspec_ttu.png'
end 
% for i=nt1:nt2
%     umax1= umax(~isnan(umax));
%     if(~isnan(umax(i)))
%         umax(i)=0.0; 
% %     end
% %       R=corrcoef(umax,umax_emp);
% %       RR=R(1,2); 
%     end
% % end
% umax1=umax; 
% for i=nt1:nt2
%     cff=isnan(umax(i))
%     cff1=int8(cff)  ;
%     umax1(i)=umax1(i)*(1-cff1)  
% end 

%figure(9)
%scatter(umax,umax_emp,'.')
% %if(~isnan((umax)))
%  R=corrcoef(umax,umax_emp);
%  RR=R(1,2)
%end 
%title(['Ucrest, correlation coeff is',num2str(RR)])
%xlabel('Direct-ADV')
%ylabel('Empirical-WH')
%print -dpng '../pngfiles/waveform_compare/waveformscatter_wh_corrceff_ucrest.png'
 
