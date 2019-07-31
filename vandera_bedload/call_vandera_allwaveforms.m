clear all ; close all ; clc; 
% Enter sediment information in meters
d50 = 0.4e-3 ; %0.2e-3;
d90 = 1.3*d50;

% Near bottom current data
% This is constant
deg2rad=pi/180.0; 
umag_curr=0.0; %.2814;% abs(0.0);
phi_curwave=0.0;% 79.92*deg2rad ;% 0.0*deg2rad;
% Zref is a reference height for computing current friction factor 
Zref=0.04 ;
% delta is the reference height at which current velocity is computed (Wave boundary layer thickness) 
delta=0.2;


% check Steve's wave form
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/9917adv_wfr.mat')

load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/9917adv_waveforms.mat'); 
% Use the detph from adv generated data 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Hrmsu','depth')
h=depth; 

%wfr=[b(1).wf(1)]

it=0 ; 
s2=length(struct(b)); % length of the wavebursts saved here
for j=1:s2
 s1=length(struct(b(j).wf)); % find the lenght of all the waveforms 
                             % stored in the structure of 1 wave
 %within each wave burst
 for i=1:s1
   wfr(i)=[b(j).wf(i)]; % check the first waveburst has 107 elements  
   umax(i+it)=[wfr(i).umax] ; 
   umin(i+it)=[wfr(i).umin] ;
   T_c(i+it)=[wfr(i).Tc]   ;
   T_t(i+it)=[wfr(i).Tt]   ;
   T_cu(i+it)=[wfr(i).Tcu] ; 
   T_tu(i+it)=[wfr(i).Ttu] ; 
   T(i+it)=[wfr(i).T] ; 
    
   h_1d(i+it)=h(j);    % Convert the dimensions of depth to what is needed ! 
   uhat(i+it)=0.5*(umax(i+it)-umin(i+it)) ;
  end
  it=s1+it ;  
end

 

% Question- Would need to use peak wave period and not bottom wave period
waveavgd_stress_term=1; 
surface_wave=1;  
current_timeperiod=0; 


for i=1:length(uhat)
   if(~isnan(uhat(i))) 
   Hs(i)=0.0; % This is a redundant input( these are hardwired to be zero because they are not even needed for direct waveforms
   % as Steve already got all the inpputs.) ..........
   
   
   [bedldx_allwaveform(i), bedldy(i), bedldtx(i), Ur(i), R(i), Beta(i)]=...........
                  vandera_function_directwaveform(i, Hs(i), T(i), h_1d(i), d50, d90, .....
                                                  umag_curr, phi_curwave, uhat(i), .....
 						                          umax(i), umin(i), ........
                                                  T_c(i), T_t(i), T_cu(i), T_tu(i), ............
 						                          Zref, delta, waveavgd_stress_term, ......
                                                  current_timeperiod, surface_wave); 
 emnd 
end  
save('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_allwaveforms.mat','bedldx_allwaveform','R','Beta','Ur')
  