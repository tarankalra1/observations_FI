clear all ; close all ; clc; 
% Use the workhorse data along with ubspecdat code with vspec
% to giev teh vandera bedload 
% Code written by Tarandeep S. Kalra and Chris Sherwood

% Enter sediment information in meters
d50 = 0.4e-3 ; %0.2e-3;
d90 = 1.3*d50;

% Near bottom current data
% This is constant
deg2rad=pi/180.0; 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/unet_9917adv_phi.mat',.......
                                               'phid','wdir','cdir','unet')

umag_curr=unet  ; % unet;  %.2814;% abs(0.0);
phi_curwave=phid;% 79.92*deg2rad ;% 0.0*deg2rad;

% Zref is a reference height for computing current friction factor 
Zref=0.1;
% delta is the reference height at which current velocity is computed (Wave boundary layer thickness) 
delta=0.1;

nt1=1; nt2=2044;
%  load('skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%             'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube')

% LOAD the observational data from workhorse from Fire Island 
url='/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'
netcdf_load(url)
Hs(:)=squeeze(wh_4061(1,1,:));
%Td(:)=squeeze(wp_peak(1,1,:));
h(:)=squeeze(hght_18(1,1,:)); % extract depth; 
initial_sensor_height = 2.11;
% Add the height to depth
h(:)=h(:)+initial_sensor_height ;
 
% VSPEC 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vspec_uhat_tr.mat','uhat_wh','Tr_wh')


% Question- Would need to use peak wave period and not bottom wave period
waveavgd_stress_term=1; 
surface_wave=1;
current_timeperiod_mod=0; 


% NOTE THAT I AM USING TBOT AND NOT PEAK WAVE PERIOD 
% 
%unet(150)
 for i=nt1:nt2
  % if(unet(i)>0) 
   % if(~isnan((i)))
  % if(~isnan(unet(i)))
    if(~isnan(Tr_wh(i)) & ~isnan(unet(i)))   
       [bedldx_wh_unet(i), bedldy_wh_unet(i), Ur(i), R(i), Beta(i), Tc_ini(i),Tc(i),Tt_ini(i),Tt(i)]=........
                              vandera_function_workhorse(i, Hs(i), Tr_wh(i), h(i), d50, d90, .....
                                                  umag_curr(i), phi_curwave(i), uhat_wh(i), .....
                                                  Zref, delta, waveavgd_stress_term, ......
                                                   surface_wave, current_timeperiod_mod) ;
                                   
    else
     bedldx_wh_unet(i)=0;
     bedldy_wh_unet(i)=0; Ur(i)=0; R(i)=0; Beta(i)=0; 
    end 
 end  
  
%   
if(current_timeperiod_mod==1);
 save('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet.mat',.....
   'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc_ini','Tc','Tt_ini','Tt')
end 

if(current_timeperiod_mod==0); 
 save('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet_timezero.mat',.....
     'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur','Tc_ini','Tc','Tt_ini','Tt')
end 

%  
% figure(1)
% subplot(2,1,1)
% plot((bedldx_wh_unet))
% subplot(2,1,2)
% plot((bedldy_wh_unet))
%  %  
% figure(2)
% plot(cumtrapz(bedldx_wh_unet))
% %  
% plot(Tc_ini)
% hold on
% plot(Tc,'k--')
