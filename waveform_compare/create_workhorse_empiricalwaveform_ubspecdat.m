clear all ; close all ; clc; 
% This code calculates the workhorse based empirical inputs for vandera formula
% It used the wave spectra (ubspecdat.m) for doing the calculations 

% LOAD the observational data from workhorse from Fire Island 
url='/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'
netcdf_load(url)
Hs(:)=squeeze(wh_4061(1,1,:));
%Td(:)=squeeze(wp_peak(1,1,:));
 h(:)=squeeze(hght_18(1,1,:)); % extract depth; 
nt1=1; nt2=2044; 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/vspec_uhat_tr.mat','uhat_wh','Tr_wh')
uhat_emp=uhat_wh;
Tbr=Tr_wh; 

%Read in the spectra file.mat for uhat and Tbr
% take the wave hegiht directly from workhorse , may be that calculation can be improved
for i=nt1:nt2
%    % [uhat_emp(i),Tbr(i)]=ubspecdat(squeeze(h_1d(i)),0.001*s(:,i)',f(:,i)',df);
      h_1d(i)=h(i); 
      % Used the representate time period 
      [r(i), phi(i), Ur_emp(i)]=skewness_params(Hs(i),Tbr(i),h_1d(i));
%      
      [Tc_emp(i), Tt_emp(i), Tcu_emp(i), Ttu_emp(i), umax_emp(i), umin_emp(i), RR_emp(i),....
          beta_emp(i)]=abreu_points(r(i), phi(i), uhat_emp(i), Tbr(i));
end 
%  end
% % 
Uw_emp=uhat_emp; 
% 
save('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec.mat',.....
     'Ur_emp','Hs','Tbr','h',.......
      'umax_emp','umin_emp','Tc_emp','Tt_emp',........
      'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp'); 
 
 %
 function [r, phi, Ur] = skewness_params(H_s, T, depth);
 %
 % Ruessink et al. provides equations for calculating skewness parameters
 % Uses Malarkey and Davies equations to get "bb" and "r"
 % Given input of H_s, T and depth
 % r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
 % phi   - skewness/asymmetry parameter                         [value]
 % Su     - umax/(umax-umin)                                    [value]
 % Au   - amax/(amax-amin)                                      [value]
 % alpha - tmax/pi                                              [value]
 
 p1=0.0;
 p2=0.857;
 p3=-0.471;
 p4=0.297;
 p5=0.815;
 p6=0.672;
 dtr = pi/180.;
 %
 % Ruessink et al., 2012, Coastal Engineering 65:56-63.
 %
 % k is the local wave number computed with linear wave theory.
 %
 k=kh_calc(T,depth)./depth     ; 
 %
 % H_s=sqrt(2.0)*H_rms
 a_w=0.5*H_s ;
 Ur=0.75*a_w*k/((k*depth)^3.0);
 %
 % Ruessink et al., 2012 Equation 9.
 %
 cff=exp( (p3-log10(Ur)) /p4);
 B=p1+((p2-p1)/(1.0+cff));
 psi=-90.0*dtr*(1.0-tanh(p5/Ur^p6));
 %
 % Markaley and Davies, equation provides bb which is "b" in paper
 % Check from where CRS found these equations
 %
 bb=sqrt(2.0)*B/(sqrt(2.0*B^2.0+9.0));
 r=2.0*bb/(bb^2.0+1.0);
 %
 % Ruessink et al., 2012 under Equation 12.
 %
 phi=-psi-0.5*pi;
 %
 %
 % Where are these asymmetry Su, Au utilized
 % recreate the asymetry
 %
 Su=B*cos(psi);
 Au=B*sin(psi);
 
 return
 end % function skewness_params
 
 %% abreu_points
 function [DTc, DTt, DTcu, DTtu, umax, umin, RR, beta] = ...
     abreu_points(r, phi, Uw, T );
 %
 %  Calculate umax, umin, and phases of asymmetrical wave orbital velocity
 %
 %  Use the asymmetry parameters from Ruessink et al, 2012
 %  to get the umax, umin and phases of asymettrical wave
 %  orbital velocity to be used by Van Der A.
 %  T_c is duration of crest
 %  T_cu Duration of accerating flow within crest half cycle
 %
 omega=2.0*pi/T;
 %
 phi_new=-phi;
 
 % Malarkey and Davies (Under equation 16b)
 P=sqrt(1.0-r*r);
 %
 % Malarkey and Davies (Under equation 16b)
 %
 b=r/(1.0+P);
 %
 % Appendix E of Malarkey and Davies
 %
 c=b*sin(phi_new);
 %
 cff1=4.0*c*(b*b-c*c)+(1.0-b*b)*(1.0+b*b-2.0*c*c) ;
 cff2=(1.0+b*b)^2.0-4.0*c*c ;
 ratio=cff1/cff2;
 if(ratio>1.0)
   ratio=1.0  ;
 end
 if(ratio<-1.0)
   ratio=-1.0 ;
 end
 tmc=asin(ratio) ;
 
 %
 cff1=4.0*c*(b*b-c*c)-(1.0-b*b)*(1.0+b*b-2.0*c*c);
 cff2=(1.0+b*b)^2.0-4.0*c*c;
 ratio=cff1/cff2;
 if(ratio>1.0)
   ratio=1.0  ;
 end
 if(ratio<-1.0)
   ratio=-1.0 ;
 end
 tmt=asin(ratio) ;
 
 
 if(tmt<0.0)
     tmt=tmt+2.0*pi;
 end
 if(tmc<0.0)
     tmc=tmc+2.0*pi;
 end
 %
 % Non dimensional umax and umin, under E5 in Malarkey and Davies
 %
 umax=1.0+c;
 umin=umax-2.0;
 %
 %       Dimensionalize
 %
 umax=umax*Uw;
 umin=umin*Uw;
 %
 % phase of zero upcrossing and downcrossing (radians)
 %
  tzu=asin(b*sin(phi_new)) ;
 tzd=2.0*acos(c)+tzu;
 %
 % MD, equation 17
 %
 RR=0.5*(1.0+b*sin(phi_new));
 %
 % MD, under equation 18
 %
  
 if(r<=0.5)
     F0=1.0-0.27*(2.0*r)^(2.1);
 else
     F0=0.59+0.14*(2.0*r)^(-6.2) ;
 end
 %
 % MD, Equation 15a,b
 %
 if(r >= 0.0 && r<0.5)
     betar_0=0.5*(1.0+r);
 elseif(r>0.5 && r<1.0)
     cff1=4.0*r*(1.0+r);
     cff2=cff1+1.0;
     betar_0=cff1/cff2;
 end
 % THIS is wrong but can be fixed later;
 betar_0=0.5 ;
 %
 % MD, Equation 18
 %
 cff=sin((0.5*pi-abs(phi_new))*F0)/sin(0.5*pi*F0) ;
 beta=0.5+(betar_0-0.5)*cff;
  %
 % MD, Table 1, get asymmetry parameterization
 % using GSSO (10a,b)
 %
 cff=sqrt(2.0*(1.0+b*b)^3.0);
 Sk=3.0*sin(phi_new)/cff;
 Ak=-3.0*cos(phi_new)/cff;
 %
 % These are the dimensional fractions of wave periods needed by Van der A eqn.
 % TSK - Check source of these equations
 %
 w=1.0/omega;
 DTc=(tzd-tzu)*w;
 DTt=T-DTc;
 DTcu=(tmc-tzu)*w;
 DTtu=(tmt-tzd)*w;
 %
 T_tu=tzd*w;
 T_cu=tzu*w;
 T_c=tmc*w;
 T_t=tmt*w;
 %
 %fprintf(fid,'R, beta %f, %f\n', RR, beta);
 return
 end % function abreu_points
 
 %% kh_calc
 function kh=kh_calc(Td,depth);
  
 %
 %  Calculate wave number from Wave period and depth
 %
 % RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
 % HR Wallingford Report TR 155, February 2006
 %
 g = 9.81;
 omega=2.0*pi./Td;
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
 
 
