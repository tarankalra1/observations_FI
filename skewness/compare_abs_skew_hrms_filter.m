clear all ; close all ; clc ;

% Comapare skewness with workhorse and directly from ADV 

nt1=1; nt2= 2044; 
% 
%nt1=680 ; nt2=730 ; 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')

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
 
 % Ursell number calcualtion from Ursell number and empirical relationships
 % to get Su, Au..
 
 rp = taran_ruessink_empirical_skewness(Ur(i));
 Su_skewness_ruess(i)=rp.Su;
 Au_skewness_ruess(i)=rp.Au;  
end 
 
%  
% 
% 
% REMOVE HRMSU <0.7 
for i=nt1:nt2
  if(Hrmsu(i)<0.7)
    Su_skewness_adv(i)=0.0;
    Su_skewness_ruess(i)=0.0; 
   end 
end
% Remove Absolute(skewness) < 0.1

for i=nt1:nt2 
   if(abs(Su_skewness_adv(i)<0.1))
     Su_skewness_adv(i)=0.0;
     Su_skewness_ruess(i)=0.0; 
   else 
     Su_skewness_adv(i)=Su_skewness_adv(i);
     Su_skewness_ruess(i)=Su_skewness_ruess(i);
  end 
end
%
% 
% % ARTIFICIALLY MAKE Sskewness to positive numbers
% Su_skewness_adv_abs=abs(Su_skewness_adv);

figure(3)
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'k');
hold on 
plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'r'); 
%hold on 
%hold on 
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2)*0.0,'linewidth',2); 
legend('skewness adv-measured','skewness-ruessink','Location','Southwest')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 'keeplimits')
ylim([0.0 0.3])
print -dpng  '-r300' '../pngfiles/Su_abs_skewness_gt_0_1_Hrmsu_gt_0_7.png'   
%  
