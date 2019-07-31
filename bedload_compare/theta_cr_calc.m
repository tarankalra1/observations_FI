function theta_cr = theta_cr_calc(d50, rhos);
%
% Critical Shields parameter from Soulsby (1997).
%
rho0=1025.0;
nu = 1.36E-6;
g = 9.81;

s=rhos/rho0;
dstar=(g*(s-1)/(nu*nu))^(1.0/3.0)*d50;
cff1=0.30/(1.0+1.2*dstar);
cff2=0.055*(1.0-exp(-0.020*dstar));
theta_cr=cff1+cff2;
%
return
end % function theta_cr_calc

