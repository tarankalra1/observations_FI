clear all ; close all ; clc; 

nt1=1; nt2= 2044; 

load('matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube') 
 
load('matfiles/vandera_bedld.mat','bedldx','R','Beta','Ur')

 figure(4) 
 subplot(2,1,1)
 plot(dn(nt1:nt2), (ur_cube(nt1:nt2)),'r--'); 
 
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 legend('measured bedld')
 
 subplot(2,1,2)
plot(dn(nt1:nt2), (bedldx(nt1:nt2)),'bo-'); 
legend('Calculated bedld')
 %ylabel('ur-mean cube') 
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % '
 print -dpng 'bedload_s                                         eparate.png'  
 
 
 
 figure(5) 
 %subplot(2,1,1)
 plot(dn(nt1:nt2), (ur_cube(nt1:nt2)),'r--'); 
 hold on
 plot(dn(nt1:nt2), (bedldx(nt1:nt2))*1000,'bo-'); 

 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 
 legend('measured bedld','1000*calculated bedld')
 
% subplot(2,1,2)
%legend('Calculated bedld')
 %ylabel('ur-mean cube') 
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % '
 print -dpng 'bedload_together.png'  