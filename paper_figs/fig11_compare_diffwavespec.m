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
% different pspec 
s_pspec(:,:)=double(pspec(1,1,:,:));  
% different vspec 
s_vspec(:,:)=double(vspec(1,1,:,:));  
% different sspec 
s_sspec(:,:)=double(sspec(1,1,:,:));  

band_width=0.015625  ;
f=double(frequency(:,1));
df=band_width ;


initial_sensor_height = 2.11;

% Add the height to depth
h_1d(:)=h_1d(:)+initial_sensor_height ;

s=(s*0.001).^2; % converted to sq.m/Hz from sq.mm/Hz 

for t=nt1:nt2; 
    s_pspec(s_pspec<0)=0.0;
    s_pspec(s_pspec>1e12)=0.0;
    % 
    s_vspec(s_vspec<0)=0.0;
    s_vspec(s_vspec>1e12)=0.0;
    %
    s_sspec(s_sspec<0)=0.0;
    s_sspec(s_sspec>1e12)=0.0;
end 

% get dn ADV DATA
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','dn'); 

nt=150 ; % corresponding to Valentine's day storm 

figure(1)

plot(dn(nt1:nt2),s_pspec(:,nt)','r');
hold on
plot(dn(nt1:nt2),s_vspec(:,nt)','b');
hold on
plot(dn(nt1:nt2),s_sspec(:,nt)','g');

xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
legend('Pspec','Vspec','Sspec')
xlabel('Time (year/month/day)')
ylabel('Wave spectrum ')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
%print('-dpng','-r300','fig2_ubr.png')

