%% puv_proc_FI - Process ADV records from Fire Island ADV
% based on puv_proc_MVCO2011.m version of Feb 13, 2009
% last revised 5/1/2019
%clear all ; close all ; clc ; 

kN1 = 30*0.018e-2 % roughness to be used on GM model
kN2 = 30*0.4e-2

advbfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advb-cal.nc'); % burstfile name
advsfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advs-cal.nc'); % statistics filename
metfn='C:\Users\ssuttles\data\FireIsland\analysis\Taran\9851met.cdf'; %nearby met buoy data
scfn='C:\Users\ssuttles\data\FireIsland\analysis\Taran\9912sc-a.cdf'; %seacat logger with temp and salinity data


ncload(advsfn); % load the statistics file

zr = 0.4; % placeholder...need to check measurement elevation here
bps.instrname = 'ADV 9917';
%bps.instrname= 'ADV 9885'; 

% you can see all of the variables in the burst file with:
% ncdisp(advbfn);
%time=ncread(advbfn,'time');
%time2=ncread(advbfn,'time2') ; 

jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

gb_first = 1; % assume first burst is good
gb_last = length(dn); % there are some bad bursts
gb = (gb_first:gb_last)';

dn = dn(gb);

brange = brange(gb); % Height of ADV transducer above boundary
brangei=interp1(dn(brange<1e30),brange(brange<1e30),dn); %interpolate bad values if possible


zoff = ncreadatt(advbfn,'/','ADVProbeSamplingVolumeOffset')/100.
zr_median = nanmedian(brangei)-zoff % this is the median ADV sample elevation
zr = brangei-zoff; % time series of ADV sample locationszr

% zr has some pretty low values, indicating that measurements were being
% made 2 cmab. We need to check this and think about it.

z_init = ncreadatt(advbfn,'u_1205','initial_sensor_height')
p_z = ncreadatt(advbfn,'P_4022','initial_sensor_height');
pdelz = p_z-(z_init-zoff); % elevation diff. between velocity and pressure- NEED to include zoff to get distance between velocity and pressure measuremens
zp = zr+pdelz; % elevation of pressure measurements [m] accounting for variable brange

zpmf=medfilt(zp,15); %smooth using median filter

%% Correct pressure data for zero offset and changes in amospheric pressure

    %load met data with barometric pressure
    baro=squeeze(ncread(metfn,'BP_915')); %barometric pressure in millibars
    dnbaro=nctime2dn(metfn);
   
    pressdb=P_4023*0.01;
    
    bproff = 0.0470; %zero offset found comparing readings in air of pressure sensor and local atmospheric pressure
    
    pressdb_ac=findp1ac(pressdb,dn,baro/100,dnbaro,bproff); %correct pressure for zero offset and atmosperic pressure changes over deployment
%% need to convert corrected pressure to depth using temp & salinity data from deployed TC data

sc.sal=squeeze(ncread(scfn,'S_41'));
sc.temp=squeeze(ncread(scfn,'T_28'));
sc.dn = nctime2dn(scfn)';

rhosw = sw_dens0(sc.sal,sc.temp);
rhoswi=interp1(sc.dn,rhosw,dn);

depac = pressdb_ac.*10^4./rhoswi./9.81; %depth of sensor
    
depth_corrected = zpmf+depac; % time series of depth [meters]
depth_corrected=depth_corrected;


save mat\9917advs_depth_corrected depth_corrected dn


%%
function p1ac=findp1ac(p1db,dn,bprdb,dnbpr,bproff)
%p1ac=p1ac(p1db,dn,bprdb,dnbpr,bproff)
%m-fcn to find atmospherically corrected pressure variable (P_1ac)
%
%INPUTS:
%
%   p1db = pressure (decibars); vector, matrix, or multi-dimensional netcdf
%
%   dn = matlab datenum time record for pressure record;
%
%   bprdb = barometric pressure from local/regional source (decibars); vector
%   only!! NOTE: BE SURE TO CONVERT UNITS FOR ATM PRESS RECORD TO DECIBARS!!!
%
%   dnbpr = matlab datenum time record for barometric pressure data in
%   same time base as netcdf file (VERY IMPORTANT!!!); vector only
%
%   bproff = offset of pressure data in air to local barometric pressure
%   p1_air-baro (decibars); can find by using scripts "findatmoffset2p1_interactive.m" or
%           "findatmoffset2p1_manual.m" or use single point calibration.
%   if p1 is a gauge pressure reading (i.e. zeroed in air) then bproff will be negative.
%   SIGN of bproff is IMPORTANT!!
%
%OUTPUT:
%
%   p1ac = sub-pressure record corrected to remove offsets and changes due
%   atmospheric  pressure.
%
%S. Suttles- 21-Sep-2016



%interpolate barometric pressure data to p1 time
bprdbi=interp1(dnbpr,bprdb,dn);
acor = bproff + bprdbi;

%convert fills to NaN in P_1
flv=1e35;
p1db(p1db==flv)=NaN;

%apply correction
p1ac=p1db(:)-acor(:);

%replace NaNs with fill value
p1ac(isnan(p1ac))=flv; 

% reshape the single vector back into the 2D sizeof P_1
p1ac=reshape(p1ac,size(p1db));

%return
end