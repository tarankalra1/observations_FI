clear all ; close all ; clc ; 

%Summary of this function goes here
%   Detailed explanation goes here
%urms=1 ;
%Uw=sqrt(urms);
%Uw=0.95;
% %T=6.5 ; depth=3.5 ; H_s=1.65;
load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
             'Hrmsu','Tr','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
% 
% % 
T=double(Tr);  % representative peak period
H_s=double(Hrmsu);  % 
Uw=double(Ubr) ; 
% 
ok = find(~isnan(depth));
% %ntime=30; 
% nt1=1; nt2=2044; 

% create a wave form
%

%urms=1; 
%Uw=sqrt(urms);
% Uw=0.95;
%T=6.5 ; depth=3.5 ; H_s=1.65;

count=0 ;
for it2=1:length(H_s)
  [r,phi,omega, umin, umax]=wave_form_parameters(H_s(it2),T(it2),depth(it2),Uw(it2));

% create a wave form 
  f=sqrt(1.0-r*r) ;

  cff=r*sin(phi)/(1.0+f) ;
  time=(linspace(0,T(it2),40));
  if(it2==1)
    T_add(it2)=0;
  else 
    T_add(it2)=T(it2) ;
  end 
  for it=1:length(time)
  
%   f=sqrt(1.0-r*r) ;

     t=time(it)  ;

%   cff=r*sin(phi)./(1.0+f) ;
%    
    cff1=sin(omega.*t)+cff ;
    cff2=1.0-r.*cos(omega.*t+phi) ;
    u(it)=Uw(it2).*f.*cff1./cff2 ;  
  
    u_tot(it+count)=u(it); 
    

%    T_tot(it+count)=time(it)+T(it2); 
  end 
  count=count+40; 
%    plot(time(it)/T(1),u,'r.')
end
save('utot_art_waveform.mat','u_tot')
%plot(ur_maj_rot(680))
%hold on
%plot(u_tot(27200))
% figure(1)
% plot(time_tot/T_tot,u_tot,'r-*')
%
% These are the dimensional fractions of wave periods needed by Van der A eqn.
% TSK - Check source of these equations
%
%w=1.0/omega;
%DTc=(tzd-tzu)*w;
%DTt=T-DTc;
%DTcu=(tmc-tzu)*w;
%DTtu=(tmt-tzd)*w;
%%
%T_tu=tzd*w;
%T_cu=tzu*w;
%T_c=tmc*w;
%T_t=tmt*w;
%
