function [T_c_new, T_cu_new, T_t_new, T_tu_new]=surface_wave_mod(Td, depth, uhat,....
                                       T_c, T_cu, T_t, T_tu)
% 
% Crest period extension for horizontal particle displacement.
% Tuning factor eps_eff = 0.55 from measurements GWK Schretlen 2010         
% Equations in Appendix B of SANTOSS Report 
% Report number: SANTOSS_UT_IR3
% Date: January 2010 
% Inputs:
% T - Time period of the waves
% depth 
% Uw- orbital velocity amplitude (uhat) 
% Time periods in crest and trough
% Outputs:
% Modified time periods in crest and trough
% 
k=kh_calc(Td,depth)/depth;
c=2*pi/(k*Td) ;         % Wave speed

eps_eff=0.55; %  r therefore eps_eff=0)

delta_Tc=eps_eff*uhat/(c*pi-eps_eff*2.0*uhat);
T_c_new=max(T_c+delta_Tc,0.0);
T_cu_new=max(T_cu*T_c_new/T_c,0.0);

delta_Tt=eps_eff*uhat/(c*pi+eps_eff*2.0*uhat);
T_t_new=max(T_t-delta_Tt,0.0);
T_tu_new=max(T_tu*T_t_new/T_t,0.0) ;

return
end 

function kh = kh_calc(Td,depth);
%
%  Calculate wave number from Wave period and depth
%
% RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
% HR Wallingford Report TR 155, February 2006
%
g = 9.81;
omega=2.0*pi/Td;
x=omega^2.0*depth/g;
%
if(x<1.0);
    y=sqrt(x);
else
    y=x;
end
%
% Iteratively solving 3 times for eqn.7 of Soulsby 1997 by using
% eqns. (12a-14)
%
t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;
%
t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;

t=tanh(y);
cff=(y*t-x)/(t+y*(1.0-t*t));
y=y-cff;
kh=y;
%
return
end % function kh_calc
