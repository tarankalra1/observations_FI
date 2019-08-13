clear all ; close all ; clc; 
% written by Tarandeep S Kalra
% to calculate mpm bedload but use the averaged waveform instead of instantaneous
% u and v
% the averaged waaveform could be the ubspecdat bassed empirical waveform 
% I tried to split the crest and trough in mpm but this is not the right
% way

nt1=1; nt2= 2044; 
%      
gamma=8; 
% make sure same d50 used for the vandera formula
d50 = 0.4e-3 ; %0.2e-3;
d90 = 1.3*d50;
rho0 = 1025.0;
g = 9.81 ;            % m/s2
rhos= 2650.0;

s=rhos/rho0; 
smgd=(s-1.0)*g*d50;
osmgd=1.0/smgd; 

smgd3=(s-1)*g*d50.^3; 
theta_cr=theta_cr_calc(d50,rhos)  ; 

% calculate for each burst 
burst_time=1 ; 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_pspec_inithght.mat',.....
      'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
       'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');


for t=nt1:nt2
   omega_tc(t)=2.0*pi/Tc_emp(t); 
   omega_tt(t)=2.0*pi/Tt_emp(t); 
   
    % Dividing the bedlaod in + ve crest cycle 
   [qb_measured_umax(t)]=func_calc_mpm (umax_emp(t),omega_tc(t),......
                        uhat_emp(t),d50,burst_time,osmgd,theta_cr,.......
                        gamma,smgd3)  ;

    % Dividing the bedlaod in - ve crest cycle 
   [qb_measured_umin(t)]=func_calc_mpm (umin_emp(t),omega_tt(t),......
                        uhat_emp(t),d50,burst_time,osmgd,theta_cr,.......
                        gamma,smgd3) ;
			
    qb_measured(t)=qb_measured_umax(t)+qb_measured_umin(t); 
end
 
clf
dt=3600 ; 

bedload_measured=cumtrapz(qb_measured) 
plot(bedload_measured*dt) ;

save('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_avgwaveform.mat','bedload_measured'); 

%save('../../matfiles/mpm_only_ss.mat','Su_skewness','qb_measured','bedload_measured')
