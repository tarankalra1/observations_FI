clear all ; close all ; clc; 

nt1=1; nt2= 2044; 
% 
% 
load('../../matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar',...
            'ur_cube','ang_rot','Au_skewness')
 
 load('../../matfiles/skewness_orbital_array.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
           'Hrmsu','Ubr','Ursell','dn','jtb_rec') 
       
    
       
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
%theta_cr=0.0 ; 

smgd3=(s-1)*g*d50.^3; 
theta_cr=theta_cr_calc(d50,rhos)  ; 

% calculate for each burst 
burst_time=8400 ; 

%jt = time+time2/(3600*24*1000);
%dn=j2dn(time,time2);
% theta_cr =0.   ; % WRONG ASSUMPTION 
 for t=1:nt2
   [qb_measured(t)]=func_calc_mpm (ur_maj_rot_array(:,t),omega_br(t),......
                        Ubr(t),d50,burst_time,osmgd,theta_cr,.......
                        gamma,smgd3) ; 
 end
  
% 
% plot(qb_measured)
%  ylim([-0.05 0.15])
%  xlim([dn(nt1) dn(nt1+300)]);
%  datetick('x',2) % '
% % hold on
% plot(bedldx(1,1:301)*3000)
clf
dt=3600 ; 

bedload_measured=cumtrapz(qb_measured) 

save('../../matfiles/mpm_only_ss.mat','Su_skewness','qb_measured','bedload_measured')
