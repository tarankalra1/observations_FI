function [qb_measured]=func_calc_mpm ( ur_maj_rot_array, omega_br,......
                                  Ubr, d50, burst_time, osmgd, theta_cr,.......
                                  gamma, smgd3); 
%calc mpm  
qb_measured=0.0 ;  
for t=1:burst_time

    if(Ubr==0)
      cff=0 ;
      cff1=0;
      fw=0.0;
    else 
     cff=omega_br./Ubr ;   %equation 8 Shawn paper
     cff1=(2.5*d50*cff).^0.2; 
     fw=exp(5.5*cff1-6.3); 
    end 
  
   theta(t)=0.5*fw.*ur_maj_rot_array(t).*abs(ur_maj_rot_array(t) )*osmgd ; 
   
   cff1(t)=(max((abs(theta(t))-theta_cr), 0.0)).^1.5;
   
% multiply with the sign of theta(t)
   cff1(t)=cff1(t)*sign(theta(t));
% %   
 %  qb_measured=gamma*cff1(t)*sqrt(smgd3)+qb_measured; 
qb_measured=gamma*cff1(t)*sqrt(smgd3)+qb_measured; 
end
qb_measured=qb_measured/8400; 
