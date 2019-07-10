clear all ; close all ; clc; 

% Cummulative skewness

clear all ; close all ; clc ;

% Comapare skewness with workhorse and directly from ADV 

nt1=1; nt2= 2044; 
% 
%nt1=680 ; nt2=730 ; 

load('../../matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness')
 Su_skewness_old=0.0;     
 for t=nt1:nt2
   Su_skewness_adv(t)=Su_skewness(t); %+Su_skewness_old; 
   Au_skewness_adv(t)=Au_skewness(t); %+Su_skewness_old; 
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
 Au_skewness_ruess(i)=rp.Au;  
end 
 
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


% ARTIFICIALLY MAKE Sskewness to positive numbers
Su_skewness_adv_abs=abs(Su_skewness_adv);
%  

meas_skew=cumtrapz(Su_skewness_adv_abs); 
calc_skew=cumtrapz(Su_skewness_ruess)  ; 

figure(1)
%subplot(2,1,1)
plot(dn(nt1:nt2),meas_skew)
%title('Cumulative measured skewness')
%xlim([dn(nt1) dn(nt2)]);
%datetick('x',2) % 
% 
% % % integrate in time 
% subplot(2,1,2)
 hold on 
plot(dn(nt1:nt2),calc_skew)
% title('Cumulative calculated skewness')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) %
 legend('Cumulative measured skewness','Cumulative calculated skewness')
 print -dpng '../../pngfiles/cumulative_skewness.png'


 