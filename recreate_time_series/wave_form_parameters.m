function [r, phi, omega, umin, umax]=wave_form_parameters(H_s,T,depth,Uw)
%
% 1. GET URSELL NUMBER 
% H_s=sqrt(2.0)*H_rms
a_w=0.5*H_s;

kh=qkhfs( 2*pi/T, depth);
k=kh/depth; 
Ur=0.75*a_w*k/((k*depth)^3.0)  ; 
%Ursell(n) = 0.75*0.5*PUV(n).Hrmsu*k(n)./(kh.^3); % RRvR Eqn. 6.


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
r=2.0*bb/(bb^2.0+1.0) ;
%
% Ruessink et al., 2012 under Equation 12.
%
phi=-psi-0.5*pi;
%
omega=2.0*pi/T;

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

