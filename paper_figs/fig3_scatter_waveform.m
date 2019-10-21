clear all ; close all ; clc; 
%
% WORKHORSE
%
% Used ubspecdat. m to obtain empirical waveform with vspec and compare with Steve's manually
% obtained waveforms.
   
 load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec_inithght.mat',.....
      'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
       'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');  
% workhorse_emp_waveform_ubspecdat_vspec_inithght.mat

load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 

%load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_wfr.mat')
load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_wfr.mat')

x_str=linspace(0,10,100); 
x_str1=linspace(-1,1,200); 
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
end


isave=1;
if(isave==1)

figure(1)
subplot(3,2,1)
scatter(umax,umax_emp,'r.')
hold on 
plot(x_str1,x_str1,'k')
%axis image
xlim([0.0 0.9])
ylim([0.0 0.9])
R=corrcoef(umax,umax_emp,'row','complete');
RR=R(1,2)
text(0.13,0.73,['r^{2} = ' (num2str(RR)) ] )

%title(['ucrest, correlation coeff is ',num2str(RR)])
%title('Ucrest')
xlim([0 0.8]) ; ylim([0 0.8])
xlabel('Direct, u_{c} (m/s)')
ylabel('Parameterized, u_{c} (m/s)')

% print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_ucrest_initsens.png'
 
 %figure(3)
subplot(3,2,2)
scatter(umin,umin_emp,'r.')
hold on 
plot(x_str1,x_str1,'k')
%axis image
xlim([-0.8 0.0])
ylim([-0.8 0.0])
R=corrcoef(umin,umin_emp,'row','complete');
RR=R(1,2)
xlim([-0.8 0.0]) ; ylim([-0.8 0.0])
%title(['utrough, correlation coeff is ',num2str(RR)])
xlabel('Direct, u_{t} (m/s)')
ylabel('Parameterized, u_{t} (m/s)')
text(-0.73,-0.1,['r^{2} = ' (num2str(RR)) ] )

% print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_utrough_initsens.png'
% 
% figure(4)
% scatter(uhat,uhat_emp,'.')
% hold on 
% plot(x_str,x_str,'k')
% axis image
% xlim([0.0 0.8])
% ylim([0.0 0.8])
% R=corrcoef(uhat,uhat_emp,'row','complete');
% RR=R(1,2)
% title(['uhat, correlation coeff is ',num2str(RR)])
% xlabel('Direct')
% ylabel('Parameterized')
%print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_uhat_initsens.png'

%figure(5)
subplot(3,2,3)
scatter(T_c,Tc_emp,'r.')
hold on 
plot(x_str,x_str,'k')
%axis image
%xlim([1.0 8.0])
%ylim([1.0 8.0])
R=corrcoef(T_c,Tc_emp,'row','complete');
RR=R(1,2)
xlim([2.0 7.0]); ylim([2.0 7.0])
%title(['Tcrest, correlation coeff is ',num2str(RR)])
xlabel('Direct, T_{c} (s)')
ylabel('Parameterized, T_{c} (s)')
text(2.8,6.8,['r^{2} = ' (num2str(RR)) ] )

%print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_tcrest_initsens.png'

%figure(6)
%subplot(4,1,4)
subplot(3,2,4)
scatter(T_t,Tt_emp,'r.')
hold on 
plot(x_str,x_str,'k')
%axis image
%xlim([1.0 8.0])
%ylim([1.0 8.0])
R=corrcoef(T_t,Tt_emp,'row','complete');
RR=R(1,2)
xlim([2.0 7.0]) ; ylim([2.0 7.0])
%title(['Ttrough, correlation coeff is ',num2str(RR)])
xlabel('Direct, T_{t} (s)')
ylabel('Parameterized, T_{t} (s)')
text(2.8,6.8,['r^{2} = ' (num2str(RR)) ] )

%print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_trough_initsens.png'

%figure(7)
%subplot(4,2,1)
subplot(3,2,5)
scatter(T_cu,Tcu_emp,'r.')
hold on 
plot(x_str,x_str,'k')
%axis image
%xlim([1.0 4.0])
%ylim([1.0 4.0])
R=corrcoef(T_cu,Tcu_emp,'row','complete');
RR=R(1,2);
xlim([1.0 4.0]) ; ylim([1.0 4.0])
%title(['Tcu, correlation coeff is ',num2str(RR)])
xlabel('Direct, T_{cu} (s)')
ylabel('Parameterized, T_{cu} (s)')
text(1.5,3.8,['r^{2} = ' (num2str(RR)) ] )

%axis equal
%print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_tcu_initsens.png'

%figure(8)
subplot(3,2,6)
scatter(T_tu,Ttu_emp,'r.')
hold on 
plot(x_str,x_str,'k')
%axis image
%scatter(T_tu,Ttu_emp,'.')
R=corrcoef(T_tu,Ttu_emp,'row','complete');
RR=R(1,2);
xlim([1.0 4.0]) ; ylim([1.0 4.0])
%title(['Ttu, correlation coeff is ',num2str(RR)])
xlabel('Direct, T_{tu} (s)')
ylabel('Parameterized, T_{tu} (s)')
text(1.5,3.8,['r^{2} = ' (num2str(RR)) ] )

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8])
%print('-dpng','-r300','fig3_scatter_waveasym.png')
%print -dpng '../pngfiles/waveform_compare/waveformscatter_whpspec_ttu_initsens.png'
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
 
