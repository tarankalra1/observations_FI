function [Tc, Tt] =current_timeperiod(unet, Ang, ...
                                  umax, umin, RR_cmg, DTc_cmg, DTt_cmg, T)
 %
 % Modify the crest and trough time periods based on current velocity
 % This function was developed by Chris Sherwood and Tarandeep Kalra
 %
 % The basis for these formulations are formed from Appendix A.3 
 % in SANTOSS report. 
 % Report number: SANTOSS_UT_IR3
 % Date: January 2010
 %  
unet_xdir=unet*cos(Ang) ;
if(RR_cmg==0.5)
  Tc = 0.5*T;
  Tt = 0.5*T;
  if(unet_xdir>=umax);
    Tc = 1*T;
    Tt = 0;
  elseif(unet_xdir<=umin);
    Tc = 0;
    Tt = 1*T;
  elseif(unet_xdir< 0 & unet_xdir > umin);
    delt = asin(-unet/umin)/pi;
    Tt = Tt*(1-2*delt);
    Tc = T-Tt;
  elseif(unet_xdir > 0 & unet_xdir < umax);
    delt = asin(unet_xdir/-umax)/pi
    Tc = Tc*(1-2*delt);
    Tt = T-Tc;
   elseif(unet_xdir==0.0)
    Tc=Tc;
    Tt=Tt; 
  end
elseif(RR_cmg>0.5)
  Tc = DTc_cmg;
  Tt = DTt_cmg;  
  if(unet_xdir>=umax);
    Tc = 1*T; 
    Tt = 0;
  elseif(unet_xdir<=umin);
    Tc = 0;
    Tt = 1*T;
  elseif(unet_xdir < 0 & unet_xdir > umin);
    delt = asin(-unet_xdir/umin)/pi;
    Tt = Tt*(1-2*delt);
    Tc = T-Tt;
  elseif(unet_xdir > 0 & unet_xdir < umax);
    delt = asin(unet_xdir/-umax)/pi ;
    Tc = Tc*(1-2*delt) 
    Tt = T-Tc;
  elseif(unet_xdir==0.0);
    Tc=Tc; 
    Tt=Tt;
  end 
end
