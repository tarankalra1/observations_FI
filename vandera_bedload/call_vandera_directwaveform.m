clear all ; close all ; clc; 
% This is based on the prototype code intended for ROMS, as of
% April 15, 2019
% code to call vandera bedload routines for a time-series 
% based on workhorse and ADV data
% Code written by Tarandeep S. Kalra and Chris Sherwood

% Enter sediment information in meters
d50 = 0.4e-3 ; %0.2e-3;
d90 = 1.3*d50;

% Near bottom current data
% This is constant
deg2rad=pi/180.0; 
umag_curr=0.0; %.2814;% abs(0.0);
phi_curwave=0.0;% 79.92*deg2rad ;% 0.0*deg2rad;
% Zref is a reference height for computing current friction factor 
Zref=0.04 ;
% delta is the reference height at which current velocity is computed (Wave boundary layer thickness) 
delta=0.2;

nt1=1; nt2= 2044;

% LOAD the observational data from workhorse from Fire Island 
url='/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'
netcdf_load(url)
Hs(:)=squeeze(wh_4061(1,1,:));
Td(:)=squeeze(wp_peak(1,1,:));
h(:)=squeeze(hght_18(1,1,:)); % extract depth; 

 % h=depth; % Depth from ADV  
% %ntime=end 
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
     [uhat(i),Tbav(i)]=ubspecfun(Hs(i),Td(i),h(i) ); 
 end
% 

% Question- Would need to use peak wave period and not bottom wave period
waveavgd_stress_term=1; 
surface_wave=1;  

% Read in Steve's waveform to get umax, umin, T_c, T_t....
%

 for i=nt1:nt2
   [bedldx(i), bedldy(i), Ur(i), R(i), Beta(i)]=vandera_function_directwaveform(i, Hs(i), Td(i), h(i), d50, d90, .....
                                                  umag_curr, phi_curwave, uhat(i), .....
						 umax(i), umin(i), ........
                                                 T_c(i), T_t(i), T_cu(i), T_tu(i), ............
						  Zref, delta, waveavgd_stress_term, surface_wave);
 
end 
 save('vandera_bedld_directwaveform.mat','bedldx','R','Beta','Ur')
