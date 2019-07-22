clear all ; close all ; clc ;

% Comapare skewness with workhorse and directly from ADV 

nt1=2000; nt2= 2028; 
% 
%nt1=680 ; nt2=730 ; 

load('../../skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
             'Hrmsu','Ubr','Ur','dn','jtb_rec')
 Su_skewness_old=0.0;     
 for t=nt1:nt2
   Su_skewness_adv(t)=Su_skewness(t); %+Su_skewness_old; 
 end 

% WORKHORSE DATA from linear wave theory
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)
 Hs(:)=squeeze(wh_4061(1,1,:));    % extract significant wave  height 
 Td(:)=squeeze(wp_peak(1,1,:));    % extract peak wave period 
 depth(:)=squeeze(hght_18(1,1,:)); % extract depth; 
 Dwave(:)=squeeze(wvdir(1,1,:)); 
 
jt = time+time2/(3600*24*1000);
dn_wh = j2dn(time,time2);
dnsb_wh = datestr(datenum((dn_wh)));
  
% h=depth; % Depth from ADV  
% %ntime=end 
for i=nt1:nt2
  if (Hs(i)>100);
    Hs(i)=0.0;
  end
  if (Td(i)>30); 
    Td(i)=0.0;
  end 

%       % CALCULATE SKEWNESS FROM SURFACE WAVES 
 omega=2.0*pi/Td(i);
 k=qkhfs(omega,depth(i))/depth(i);
 a_w=0.5*Hs(i);
 Ur(i)=0.75*a_w*k/((k*depth(i))^3.0);   
 
 rp = taran_ruessink_asym(Ur(i));
 Su_skewness_ruess(i)=rp.Su;
end
%  
% figure(2)
% plot(Ur(nt1:nt2),Su_skewness_ruess(nt1:nt2),'.')
% % hold on 
% % plot(Ur(nt1:nt2),Hrmsu(nt1:nt2),'.')
% % 
% xlabel('Ursell number')
% ylabel('Su-ruessink') 
%  
 
 figure(3)
 plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'linewidth',2);
 hold on 
 plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'linewidth',2); 
 legend('skewness adv-measured','skewness-ruessink')
 xlim([dn(nt1) dn(nt2)]);
 xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % '
% print -dpng 'Su_skewness_negative_compare.png'  
 
%   datetick('x',2) % 'keeplimits')
% %  legend('adv-Ubr from burst data','workhorse-linear wave theory','workhorse-spectra')
% % ylabel('Ubr')
% 
% REMOVE HRMSU <0.5 
for i=nt1:nt2
  if(Hrmsu(i)<0.5)
    Su_skewness_adv(i)=0.0;
    Su_skewness_ruess(i)=0.0; 
    Hrmsu(i)=0.0;
  else 
    %Su_skewness_adv(i)=0.0;
    %Su_skewness_ruess(i)=0.0; 
    %Hrmsu(i)=0.0;
  end 
end
% % 
figure(4)
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'linewidth',2);
hold on 
plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'linewidth',2); 
hold on 
%hold on 
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2)*0.0,'linewidth',2); 
legend('skewness adv-measured','skewness-ruessink')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 'keeplimits')
%print -dpng 'Su_skewness_negative.png'  

figure(5)
plot(Dwave(nt1:nt2),Su_skewness_adv(nt1:nt2),'r.','linewidth',2);
%print -dpng 'Su_skewness_negative_dwave.png'
%hold on
%plot(dn(nt1:nt2),Td(nt1:nt2)/10,'linewidth',2);

