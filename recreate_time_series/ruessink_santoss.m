% Ruessink et al. provides equations for calculating skewness parameters
% Uses Malarkey and Davies equations to get "bb" and "r"
% Given input of H_s, T and depth
% r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
% phi   - skewness/asymmetry parameter                         [value]
% Su     - umax/(umax-umin)                                    [value]
% Au   - amax/(amax-amin)                                      [value]
% alpha - tmax/pi   
clear all ; clc ; close all ;

urms=0.3 ; 
Uw=sqrt(urms);
% Uw=0.95;
T=6 ; depth=5 ; H_s=2;


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
k=kh_calc(T,depth)/depth;
%
% 1. GET URSELL NUMBER 
% H_s=sqrt(2.0)*H_rms
a_w=0.5*H_s;
Ur=0.75*a_w*k/((k*depth)^3.0);  
%
% Ruessink et al., 2012 Equation 9.
%
%2. TOTAL NON LINEARITY 
%
cff=exp( (p3-log10(Ur)) /p4);
B_old=p1+((p2-p1)/(1.0+cff));

% The equation of B is different in SANTOSS thesis than the paper 
B=p2/(1.0+cff) ; 
psi=-90.0*dtr*(1.0-tanh(p5/Ur^p6)) ;

%
% Markaley and Davies, equation provides bb which is "b" in paper
% Check from where CRS found these equations
%
% THIS EQUATION IS DIFFERENT 
r=0.0517*B^3 - 0.409*B^2 + 1.0853*B - 0.0099;


%bb=sqrt(2.0)*B/(sqrt(2.0*B^2.0+9.0));
%r_old=2.0*bb/(bb^2.0+1.0);


%
% Ruessink et al., 2012 under Equation 12.
% PHASE IS RELATED TO PHI

phi=-psi-0.5*pi   
if(-pi<phi<0.0)
  r=-r;
%elseif(-pi<phi<-0.5*pi)
%  r=-r;
%else
  r=r;
end
r=0.8; phi=-5.0*pi/12.0 ; 
%
% Where are these asymmetry Su, Au utilized
% recreate the asymetry
%
Su=B*cos(psi);
Au=B*sin(psi);

omega=2*pi/T;

% create a wave form 
f=sqrt(1.0-r*r) 
% 
time=(linspace(0,T,40)); 
cff=r*sin(phi)/(1.0+f);

for it=1:length(time)
   t=time(it)  ; 
   
   cff1=sin(omega*t)+cff; 
   cff2=1.0-r*cos(omega*t+phi);
   u(it)=Uw*f*cff1./cff2;
   
   cff1=(r*r/(1.0+f))*sin(phi)*sin(omega*t+phi);
   cff2=cos(omega*t)-r*cos(phi)-cff1;
   cff3=(1.0-r*cos(omega*t+phi)).^2 ;
   
   a(it)=Uw*omega*f*cff2./cff3; 
end
figure(1)
plot(time/T,u,'r-*')
%figure(2)
%plot(time/(2*pi),a,'b-*')

