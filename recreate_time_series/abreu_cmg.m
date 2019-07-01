clear all ; close all ; clc ; 

%Summary of this function goes here
%   Detailed explanation goes here
urms=1 ;
Uw=sqrt(urms);
%Uw=0.95;
T=6.5 ; depth=3.5 ; H_s=1.65;


kh=qkhfs( 2*pi/T, depth);
k=kh/depth; 
%k=kh_calc(T,depth)/depth;
%
% 1. GET URSELL NUMBER 
% H_s=sqrt(2.0)*H_rms
a_w=0.5*H_s;
Ur=0.75*a_w*k/((k*depth)^3.0);

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
% Where are these asymmetry Su, Au utilized
% recreate the asymetry
%
Su=B*cos(psi);
Au=B*sin(psi);

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
RR=0.5*(1.0+b*sin(phi_new)) 
%
% MD, under equation 18
%
if(r<=0.5)
    F0=1.0-0.27*(2.0*r)^(2.1);
else
    F0=0.59+0.14*(2.0*r)^(-6.2);
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
%
% MD, Equation 18
%
cff=sin((0.5*pi-abs(phi_new))*F0)/sin(0.5*pi*F0);
beta=0.5+(betar_0-0.5)*cff 
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

% create a wave form 
f=sqrt(1.0-r*r) ;

cff=r*sin(phi)/(1.0+f) 
time=(linspace(0,T,40)); 

omega
for it=1:length(time)
   t=time(it)  ; 
   
   cff1=sin(omega*t)+cff ; 
   cff2=1.0-r*cos(omega*t+phi);
   u(it)=Uw*f*cff1./cff2
   
   cff1=(r*r/(1.0+f))*sin(phi)*sin(omega*t+phi);
   cff2=cos(omega*t)-r*cos(phi)-cff1;
   cff3=(1.0-r*cos(omega*t+phi)).^2 ;
   
   a(it)=Uw*omega*f*cff2./cff3; 
end
figure(1)
plot(time/T,u,'r-*')