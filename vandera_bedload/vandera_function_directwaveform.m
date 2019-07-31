% Function that takes ucrest, utrough, Tcrest, Trough, Tcu, Ttu
% Written by T.S. Kalra, July 19, 2019

function [bedld_x, bedld_y, bedld_tx, Ur, RR, beta]=vandera_function_directwaveform(time, Hs, Td, depth, d50, d90, .....
                                                 umag_curr, phi_curwave, uhat, umax, umin, .....
						                                T_c, T_t, T_cu, T_tu, ............
                                                 Zref, delta, waveavgd_stress_term, ......
                                                 current_timeperiod, surface_wave);					                      
%Summary of this function goes here
%   Detailed explanation goes here
eps_eff=1.0; 
rho0 = 1025.0;
g = 9.81 ;            % m/s2
vonKar = 0.41 ;       % non-dimensional
nu = 1.36E-6  ;       % kinematic viscosity m2 s-1
% pi = 3.14159265358979323846_r8 % use Matlab constant
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;
rhos= 2650.0;
smgd=(rhos/rho0-1.0)*g*d50;
osmgd=1.0/smgd;

% uhat is directly coming from data unlike COAWST uhat=urms*sqrt(2.0);
ahat=uhat*Td/(2.0*pi);
k=kh_calc(Td,depth)/depth;     % Wave number 
c_w=2*pi/(k*Td) ;              % Wave speed

%
% VA-2013 equation 1 is solved in 3 sub-steps
%-----------------------------------------------------------------------
umax=umax*uhat; % Uw is the wave orbital velocity  
umin=umin*uhat; 

uhat_c=umax; 
uhat_t=umin; 

uc_r=0.5*sqrt(2.0)*uhat_c;
ut_r=0.5*sqrt(2.0)*uhat_t;
%
RR=umax/(umax-umin)      ; 

a_w=0.5*Hs ;
Ur=0.75*a_w*k/((k*depth)^3.0);

beta=0.0;  %THIS IS AN ASSUMPTION BECUASE OF DEEP WATER 

[T_c, T_t, T_cu, T_tu, eta, udelta, alpha, ksw, tau_wRe]=..... 
                         full_wave_cycle_stress_factors(d50, d90, osmgd, .....
                             Td, depth, T_c, T_t, T_cu, T_tu, RR, ......
                     umag_curr, phi_curwave, Zref, delta, umax, umin, uhat, ahat);   
%
%-----------------------------------------------------------------------
% 2. Bed shear stress (Shields parameter) for Crest half cycle
%    alpha VA2013 Eqn. 19
%-----------------------------------------------------------------------
%
wavecycle=1.0;                     
[dsf_c, theta_cx, theta_cy, mag_theta_c]=half_wave_cycle_stress_factors(T_cu, T_c,......
                   uc_r, uhat_c, udelta, phi_curwave, ..........
                   alpha, fd, ahat, ksw, tau_wRe);
%
%-----------------------------------------------------------------------
% 2. Bed shear stress (Shields parameter) for Trough half cycle
%    alpha VA2013 Eqn. 19
%-----------------------------------------------------------------------
%
wavecycle=-1.0;
[dsf_t, theta_tx, theta_ty, mag_theta_t]=half_wave_cycle_stress_factors(T_tu, T_t,......
                           ut_r, uhat_t, udelta, phi_curwave,....
                           alpha, fd, ahat, ksw, tau_wRe);
%%
%-----------------------------------------------------------------------
% 3. Compute sediment load entrained during each crest half cycle
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
%      Crest half cycle
%-----------------------------------------------------------------------
%
wavecycle=1.0;
[ om_cc, om_tc ]= sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta, dsf_c,                           ...
    T_c, T_cu, uhat_c, mag_theta_c);
%
%-----------------------------------------------------------------------
%       Trough half cycle
%-----------------------------------------------------------------------
%
wavecycle=-1.0;
[om_tt, om_ct] = sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta, dsf_t,                           ...
    T_t, T_tu, uhat_t, mag_theta_t);
%-----------------------------------------------------------------------
% VA2013  Use the velocity-load equation 1.
% Non-dimensional net transport rate
%-----------------------------------------------------------------------
%
%
cff1=0.5*T_c/(T_cu);
cff2=sqrt(mag_theta_c)*(om_cc+cff1*om_tc);;
%
cff3=theta_cx/mag_theta_c;
bedld_cx=cff2*cff3 ;
%
cff3=theta_cy/mag_theta_c;
bedld_cy=cff2*cff3 ;
% 
cff1=0.5*T_t/(T_tu);
cff2=sqrt(mag_theta_t)*(om_tt+cff1*om_ct);
 
cff3=theta_tx/mag_theta_t;
bedld_tx=cff2*cff3 ;
%
cff3=theta_ty/mag_theta_t;
bedld_ty=cff2*cff3 ;
%
% The units of these are m2 sec-1
% bed_frac, rhos multiplied
%
smgd_3=sqrt((rhos/rho0-1.0)*g*d50.^3.0);

bed_frac=1.0;
bedld_x=bed_frac*(bedld_cx*T_c+bedld_tx*T_t)/Td;
bedld_y=bed_frac*(bedld_cy*T_c+bedld_ty*T_t)/Td;
 
bedld_x=bed_frac*smgd_3*(bedld_x) ;
bedld_y=bed_frac*smgd_3*(bedld_y);

%% sandload_vanderaa
function [om_ii, om_iy]=sandload_vandera(wavecycle,...
    Hs, Td,  depth, RR,                 ...
    d50, rhos, c_w,                     ...
    eta, dsf,                           ...
    T_i, T_iu, uhat_i, mag_theta_i);
%
% VA2013 Text under equation 37
%
m=11.0; n=1.2; alpha=8.2;
xi=1.7; % Based on Santoss_core.m
%
%
% Find settling velocity based on Soulsby (1997).
% VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
%
w_s=w_s_calc(0.8*d50, rhos);
 
%
% VA2013 Equation 29, for crest cycle
%  
ws_eta=w_sc_calc(Hs, Td, depth, RR, w_s, eta);
ws_dsf=w_sc_calc(Hs, Td, depth, RR, w_s, dsf);
if(wavecycle==1.0);
    w_sc_eta=max(w_s+ws_eta,0.0);
    w_sc_dsf=max(w_s+ws_dsf,0.0);
end
%
% VA2013 Equation 30, for trough cycle
%
% Not allowed it to get to zero because it gets divided 
% leading to an infinite phase lag 
%  0.36*w_s is from Shant and Dong
% 
if(wavecycle==-1.0);
    w_sc_eta=max(w_s-ws_eta,0.0) ;
    w_sc_dsf=max(w_s-ws_dsf,0.0) ; 
end
%
% VA2013 Equation 33, Phase lag parameter
%cff=(1.0-wavecycle*xi*uhat_i)/c_w;
% 
cff=(1.0-((wavecycle*xi*uhat_i)/(c_w)));
%   
cff1_eta=(1.0/(2.0*(T_i-T_iu)*w_sc_eta))  ;
cff1_dsf=(1.0/(2.0*(T_i-T_iu)*w_sc_dsf))  ; 
%
% For ripple regime 
% CRS like this:
P=alpha*dsf*cff*cff1_dsf ; 
if(eta>d50)
    P=alpha*eta*cff*cff1_eta ;
end

% TODO CRS - this is not the same eps_eff that I put into the main routine
%eps_eff=(dsf/d50)^0.25; 
eps_eff=1.0; 
theta_ieff=eps_eff*mag_theta_i ;
%
% Find critical Shields parameters based on Soulsby (1997).
%
theta_cr=theta_cr_calc(d50, rhos);
%
% Sand load entrained in the flow during each half-cycle
% 
%theta_ieff 
%theta_diff=max((theta_ieff-theta_cr),0.0) ;
% NEW CHANGE TSK
theta_diff=max((abs(theta_ieff)-theta_cr),0.0);

om_i=m*(theta_diff)^n  ;

% NEW CHANGE TSK 
om_i=om_i*sign(theta_ieff);           
%
% VA2013 Equation 23-26, Sandload entrained during half cycle
%  
if(P<=1.0)
    om_ii=om_i;
    om_iy=0.0 ;
else;
    om_ii=om_i/P; 
    cff=1.0/P;
    om_iy=om_i*(1.0-cff) ;
end
 
%if(wavecycle==1)
%elseif(wavecycle==-1)
%end 
return
end % function sandload_vandera

%% w_sc_calc
%
function worb = w_sc_calc(Hs, Td, depth, RR, w_s, zws);
%
% Second order Stokes theory to get vertical velocity of water particle
% at a given elevation based on santoss_core.m
%
worb1=pi*Hs*zws/(Td*depth);
worb2=worb1*2.0*(RR+RR-1.0);
%
%  Using the SANTOSS model formulation
%
cff=1.0/8.0;
worb=cff*worb1*sqrt(64.0-(-worb1+...
    sqrt(worb1^2.0+32.0*...
    worb2^2.0))^2.0/(worb2^2.0))+...
    worb2*sin(2.0*acos(cff*(-worb1+...
    sqrt(worb1^2.0+32.0*worb2^2.0))/worb2));

return
end % function w_sc_calc
%%
function [T_c, T_t, T_cu, T_tu, eta, udelta, alpha, ksw, tau_wRe]=....... 
    full_wave_cycle_stress_factors(d50, d90, osmgd, .....
             Td, depth, T_c, T_t, T_cu, T_tu, RR, ......
             umag_curr, phi_curwave, Zref, delta, umax, umin, uhat, ahat);
%
% Function returns
% eta-    ripple height
% udelta- current velocity at the wave boundary layer
% fd- current friction factor  
% tau_wRe- Wave averaged Reynolds stress
% T_c, T_t, T_cu, T_tu- Updated time periods in half cycles based on current velocity
% 
tol = 0.001;
total_iters=10;
%
% maximum mobility number at crest and trough
% For irregular waves, use Rayleigh distributed u(1/10) value
% VA, text under equation Appendix B.4
%
psi=(1.27*uhat)^2*osmgd;
%
% Use Appendix B eqn B.1 and B.2 to get ripple height and length
%
[eta, lambda] = ripple_dim(psi, d50); 
%
eta=eta*ahat;
lambda=lambda*ahat;

% VA2013 Eqn. 19:
%
% Initiliaze with theta_timeavg=0 and theta_hat_i=theta_timeavg
% The 
% This loop is for the full cycle 
theta_timeavg=0.0;
theta_timeavg_old=0.0;

  for iter=1:total_iters
%     % These agree 
    %   
    % Calculate wave related bed roughness from VA2013 A.5
    %   
    ksw=ksw_calc(d50, mu_calc(d50), theta_timeavg, eta, lambda);
    %   
    % Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4
    %   
    fw=fw_calc(ahat, ksw); 
    %   
    % Calculate current-related bed roughness from VA2013 Appendix A.1
    %   
    ksd=ksd_calc(d50, d90, mu_calc(d50), theta_timeavg, eta, lambda);
    %   
    % Calculate full-cycle current friction factor from VA2013 Eqn. 20
    % Within COAWST use bustr, bvstr to get delta 
    % In the standalone matlab code, delta is an input 
    %   
    % Use fd (full cycle current friction to get udelta)
    % udelta=current velocity at the top of the wave boundary layer
    %   
    fd=fd_calc(umag_curr, Zref, ksd) ;
    
    ustarc=(0.5*fd).^0.5.*umag_curr;      %friction velocity [m/s]
    %    
    udelta=max( ((ustarc/0.4)*log(30.0*delta/ksd)),0.00000001 );
    %     
    % Calculate Time-averaged absolute Shields stress VA2013 Appendix Eq. A.3
    % 
    theta_timeavg=osmgd*(0.5*fd*udelta^2.0+...
                         0.25*fw*uhat^2.0) ;
    % 
    if(abs(theta_timeavg-theta_timeavg_old) < tol);
        break
    end
    if((abs(theta_timeavg-theta_timeavg_old) >= tol)&&(iter==total_iters));
 %       fprintf(1,'Warning...stress calcs did not converge.\n');
    end
    theta_timeavg_old=theta_timeavg;
  end
  
%
% Calculate wave Reynolds stress from full cycle wave and friction factor
% that were formed from the iterative cycle, VA2013, Eqn.22
% 
alpha=udelta/(udelta+uhat);
fwd=alpha*fd+(1.0-alpha)*fw;
%
 
k=kh_calc(Td,depth)/depth;     % Wave number
c_w=2*pi/(k*Td);               % Wave speed
alpha_w=0.424;
%
if(waveavgd_stress_term==1)
 tau_wRe=rho0*fwd*alpha_w*uhat^3.0/(2.0*c_w);
else
 tau_wRe=0;
end 
%
% Compute the change in time period based on converged udelta (current velocity at 
% wave boundary layer)
%
% Calculate the time period based on udelta     
% 
if(current_timeperiod==1)
  [T_c, T_t]=current_timeperiod(udelta, phi_curwave, umax, umin, RR, T_c, T_t, Td);
end 
% Calculate the effect of surface waves 
%
if(surface_wave==1)
 [T_c, T_cu, T_t, T_tu]=surface_wave_mod(Td, depth, uhat, T_c, T_cu, .....
                                        T_t, T_tu);
end
%
return
end % function full_wave_factors

function [dsf, theta_ix, theta_iy, mag_theta_i]=half_wave_cycle_stress_factors(T_iu, T_i,......
                      ui_r, uhat_i, udelta, phi_curwave, alpha, fd, ahat, ksw, tau_wRe)
%
% For each half cycle the function returns:
% dsf - Sheet flow thickness
% theta_ix - Shields parameter in x dir. 
% theta_iy - Shields parameter in y dir. 
% mag_theta_i - Magnitude of Shields parameter
%
% Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
%
fw_i=fwi_calc(T_iu, T_i, ahat, ksw);
%
% Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
% Different for crest and trough
% 
fwd_i=alpha*fd+(1.0-alpha)*fw_i;
%
% Update theta_hat_i based on crest/trough amplitude uhat Eqn. C.2
% 
theta_hat_i=0.5*fwd_i*uhat_i^2*osmgd;
%
% Sheet flow thickness VA2013 Appendix C C.1
% Update from converged value of theta_hat_i 
%
dsf=dsf_calc(d50, theta_hat_i); %this dsf is in m
%
% Calculated the velocity magnitude based on representative velocities
% equation 12 from Van der A, 2013 
%
%-----------------------------------------------------------------------
% Get the representative trough half cycle water particle velocity
%    as well as full cycle orbital velocity and excursion
%-----------------------------------------------------------------------
%
ui_rx=ui_r+udelta*cos(phi_curwave);
ui_ry=udelta*sin(phi_curwave);
mag_ui=sqrt(ui_rx.*ui_rx+ui_ry.*ui_ry);
%
% VA2013-Magnitude of Shields parameter Eqn. 17
%
mag_theta_i=0.5*fwd_i*mag_ui.^2*osmgd;
%
%-----------------------------------------------------------------------
% Shields parameter VA2013 Eqn 15
% rho0 is required for non-dimensionalizing
%-----------------------------------------------------------------------
%
theta_ix=mag_theta_i*ui_rx/(mag_ui)+(tau_wRe*osmgd/rho0);
theta_iy=mag_theta_i*ui_ry/(mag_ui);
%
mag_theta_i=sqrt(theta_ix*theta_ix+theta_iy*theta_iy);

return
end % function half_wave_cycle_factor

%% ripple_dim
function [eta, lambda]=ripple_dim(psi, d50);
%
% Calculate ripple dimensions of O'Donoghue et al. 2006
% based on VA2013 Appendix B
% returns eta-Ripple length and lambda-Ripple length
d50_mm=0.001*d50;
if(d50_mm<0.22)
    m_eta=0.55;
    m_lambda=0.73;
elseif(d50_mm>0.22 && d50_mm<0.30)
    m_eta=0.55+(0.45*(d50_mm-0.22)/(0.30-0.22));
    m_lambda=0.73+(0.27*(d50_mm-0.22)/(0.30-0.22));
else
    m_eta=1.0;
    m_lambda=1.0;
end
%
% Smooth transition between ripple regime and bed sheet flow regime
%
if(psi<=190.0)
    n_eta=1.0;
elseif(psi>190.0 && psi<240.0)
    n_eta=0.5*(1.0+cos(pi*(psi-190.0)/(50.0)));
elseif(psi >= 240.0)
    n_eta=0.0;
end
n_lambda=n_eta;
%
eta=max(0.0,m_eta*n_eta*(0.275-0.022*psi^0.42));
lambda=max(0.0,m_lambda*n_lambda*(1.97-0.44*psi^0.21));
%
return
end % function ripple_dim

%% theta_cr_calc        %
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


%% w_s_calc
function ws = w_s_calc(d50, rhos);
%
% Critical Shields parameter from Soulsby (1997).
% Dynamics of Marine Sands
rho0=1025.0;
nu = 1.36E-6;
g = 9.81;

s=rhos/rho0;
dstar=(g*(s-1)/(nu*nu))^(1.0/3.0)*d50;
cff=nu/d50;
cff1=10.36;
ws=cff*(sqrt(cff1*cff1+1.049*dstar^3.0)-cff1);
%
return
end % function w_s_calc


%% mu_calc
function mu = mu_calc(d50);
%
% Calculate bed roughness factor based on grain size
% VA2013 Appendix A., required for current related bed roughness
% and wave related bed roughness.
d50_mm=d50*1000 ;
%
if(d50_mm<=0.15)
    mu=6.0;
elseif(d50_mm>0.15 && d50_mm<=0.20)
    mu=6.0-5.0*((d50_mm-0.15)/(0.2-0.15)) 
elseif(d50_mm>0.20)
    mu=1.0;
end
return
end % function mu_calc


%% ksd_calc
function ksd = ksd_calc(d50, d90, mu, theta_timeavg, eta, lambda);
%
% Calculate current-related bed roughness from VA2013 Appendix A.1.
%
%eta=max(eta,d50);
eta=eta;
lambda=max(lambda,d50);
ripple_fac=0.4*eta^2.0/lambda;
ksd=max( 3.0*d90,...
    d50*(mu+6.0*(theta_timeavg-1.0)) )+...
    ripple_fac;
%
return
end % function ksd_calc

%% ksw_calc
function ksw = ksw_calc(d50, mu, theta_timeavg, eta, lambda);
%
% Calculate wave related bed roughness from VA2013 Eqn. A.5.
%
%eta=max(eta,d50);
eta=eta;
lambda=max(lambda,d50);
ripple_fac=0.4*eta^2.0/lambda;
ksw = max( d50,...
    d50*(mu+6.0*(theta_timeavg-1.0)) )...
    +ripple_fac;
%
return
end % function ksw_calc

%% fw_calc
function fw = fw_calc(ahat, ksw);
%
% Calculate full-cycle wave friction factor from VA2013 Eqn. A.4.
%
ratio=ahat/ksw;
if(ratio>1.587);
    fw = 0.00251*exp(5.21*(ratio)^(-0.19));
else
    fw = 0.3;
end
%
return
end % function fw_calc

%% fd_calc
function fd = fd_calc(umag_curr, Zref, ksd);
%
% Calculate current related friction factor VA2013 Eqn. 20
% Assuming logarithmic velocity profile.
%
  if(umag_curr==0.0)
    fd=0.0;
  else 
    von_k = 0.41;
    fd = 2.0*(von_k/log(30.0*Zref/ksd))^2.0;
  end 
return
end % function fd_calc


%% fwi_calc
function fwi = fwi_calc(T_iu, T_i, ahat, ksw);
%
% Wave friction factor for wave and crest half cycle VA2013 Eqn. 21.
%
c1=2.6;
ratio=ahat/ksw;
if(ratio>1.587)
    cff=(2.0*T_iu/T_i)^c1;
   % cff1=(cff*ratio)
    fwi=0.00251*exp(5.21*(cff*ratio)^(-0.19)) ;
else
    fwi = 0.3;
end
%
return
end % function fwi_calc


%% dsf_calc
function dsf = dsf_calc(d50, theta_i);
%
% Sheet flow thickness VA2013 Appendix C.1.
%
d50_mm=d50*1000.0;
if(d50_mm<=0.15)
    cff=25.0*theta_i;
elseif(d50_mm>0.15 & d50_mm<0.20)
    cff=25.0-(12.0*(d50_mm-0.15)/0.05);
elseif(d50_mm >= 0.20)
    cff=13.0*theta_i;
end
dsf = max(d50*cff,d50);
%
return
end % function dsf_calc
%
%  end of functions for step 2 for Shear stress formulation
%
%% kh_calc
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

%% checkvals
function checkvals(name, val1, val2, tol);
if (exist('tol')~=1);tol = 1e-4; end
if( abs(val1-val2)>tol);
   % fprintf('WARNING: Diff. values for %s: %f %f\n',name,val1,val2);
end
return
end

end

