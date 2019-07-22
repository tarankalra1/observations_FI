clear all ; close all ; clc; 
% This is based on the prototype code intended for ROMS, as of
% April 15, 2019
% code to call vandera bedload routines for a time-series 
% based on workhorse and ADV data
% Code written by Tarandeep S. Kalra and Chris Sherwood


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
% Read in Steve's waveform to get umax, umin, T_c, T_t....
%
load('/media/taran/DATADRIVE2/Obs_data/matfiles/9917adv_wfr.mat')

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
 for i=nt1:nt2
      if(~isnan(uhat(i)))  
      end
 end
 
