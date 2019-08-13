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

theta_cr=0.047 ; 
for t=nt1:nt2
   ahat=uhat_emp(t)*Tbr(t)/(2.0*pi); 
   ksw=2.5*d50 ; % Is this correct ? d50
   ratio=ahat/ksw ;
  
   if(ratio>1.587)
     fw=0.00251*exp(5.21*(ratio).^(-0.19)); 
   else
     fw=0.3                               ; 
   end 
 
   tauw=0.5*fw*uhat_emp(t).^2  ; 
   theta_w=tauw*osmgd          ;       

   cff1(t)=(max((abs(theta_w)-theta_cr), 0.0)).^1.5;

% multiply with the sign of theta(t)
%   cff1(t)=cff1(t)*sign(theta(t));
% %
   qb_measured(t)=gamma*cff1(t)*sqrt(smgd3) ; 
end
  
dt=3600 ; 

bedload_measured=cumtrapz(qb_measured) ;
plot(bedload_measured*3600) ;
% % 
%  save('/media/taran/DATADRIVE2/Obs_data/matfiles/mpm_avgwaveform.mat',.......
%                         'qb_measured','bedload_measured'); 

%save('../../matfiles/mpm_only_ss.mat','Su_skewness','qb_measured','bedload_measured')
