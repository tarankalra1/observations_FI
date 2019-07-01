clear all ; close all ; clc ;
% Get ubr from ADV and compare with the ubr from workhorse
%% LOAD ADV DATA
%load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/crs_adv_9917.mat') 
     %save('crs_adv_9885.mat','dn','jtb_rec','depth','Hrmsu','Tr','Ubr')
load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
%dnsb_adv = datestr(datenum(gregorian(jtb_rec)))  ;
dn_adv=dn; 
dnsb_adv = datestr(datenum(dn_adv))   ;

ok = find(~isnan(depth));
%ntime=30; 
nt1=1; nt2=2044; 

figure(1)
plot(dn(nt1:nt2),Ubr(nt1:nt2),'linewidth',2); 
xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % 'keeplimits')
 
%% Workhorse data that is using wave spectra to get uwave_rms
load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/ubr_from_spectra.mat','ubr','Tbr')
ubr_spectra_wh=ubr; 

%% WORKHORSE DATA to get uwave_rms 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)
 Hs(:)=squeeze(wh_4061(1,1,:));
 Td(:)=squeeze(wp_peak(1,1,:));
 h(:)=squeeze(hght_18(1,1,:)); 
 
jt = time+time2/(3600*24*1000);
dn_wh = j2dn(time,time2);
dnsb_wh = datestr(datenum((dn_wh)));
% h=depth ; 
% h=isnan(depth); 
% h=depth; % Depth from ADV  
% %ntime=end 
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
     [ubr_linear_wh(i),Tbav(i)]=ubspecfun( Hs(i),Td(i),h(i) ); 
 end
% 
 hold on 
 plot(dn_wh(nt1:nt2),ubr_linear_wh(nt1:nt2),'linewidth',2);
 hold on 
 plot(dn_wh(nt1:nt2),ubr_spectra_wh(nt1:nt2),'linewidth',2); 
 xlim([dn_wh(nt1) dn_wh(nt2)]);
  datetick('x',2) % 'keeplimits')
  legend('adv-Ubr from burst data','workhorse-linear wave theory','workhorse-spectra')
  ylabel('Ubr')
 print -dpng 'ubr_plot.png'   
 