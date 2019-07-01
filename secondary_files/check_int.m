clear all ; close all ; clc; 
url='/media/taran/DATADRIVE2/Obs_data/data_netcdf/9917advb-cal.nc' ;
netcdf_load(url);


jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);
%vel = [0 .45 1.79 4.02 7.15 11.18 16.09 21.90 29.05 29.05 ...
%29.05 29.05 29.05 22.42 17.9 17.9 17.9 17.9 14.34 11.01 ...
%8.9 6.54 2.03 0.55 0];
%time = 0:24;

%vel_int=cumtrapz(vel)

%plot(vel_int)
