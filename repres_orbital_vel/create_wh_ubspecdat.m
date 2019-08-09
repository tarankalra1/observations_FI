% written by Tarandeep S Kalra with help from Steve Suttles
% it calls the output of ADV based Hrmsu, ubr, Tr and compares
% with the workhorse data to get the work horse based wave orbital velocity and representative
% time period 

clear all ; close all ; clc; 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% convert work horse surface data of wave energy spectra to get ubr and
% Tbr..

netcdf_load(wh)
nt1=1; nt2=2044; 

Hs(:)=squeeze(wh_4061(1,1,:));
h_1d(:)=double(hght_18(1,1,:));
s(:,:)=double(vspec(1,1,:,:));  

band_width=0.015625  ;
f=double(frequency(:,1));
df=band_width ;


initial_sensor_height = 2.11;

% Add the height to depth
h(:)=h(:)+initial_sensor_height ;

s=(s*0.001).^2; % converted to sq.m/Hz from sq.mm/Hz 

for t=nt1:nt2; 
    s(s<0)=0.0;
    s(s>1e12)=0.0;
    % 
    [ubr(t),Tbr(t)]=ubspecdat(squeeze(h_1d(t)),s(:,t)',f(:,1)',df);
    Hs_wh(t)=Hs(t);   
    Hs_wh(Hs_wh>100)=0.0;
end 
uhat_wh=ubr;
Tr_wh=Tbr;

save('/media/taran/DATADRIVE2/Obs_data/matfiles/vspec_uhat_tr_inithght.mat','uhat_wh','Tr_wh')
