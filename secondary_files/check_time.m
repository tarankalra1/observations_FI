clear all ; close all ; clc; 


%load('matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%            'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar',...
%            'ur_cube','ang_rot','Au_skewness')
        
%load('matfiles/ubr_from_spectra.mat');         


load('matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
             'Hrmsu','Ubr','Ur','dn','jtb_rec')

nt1=1; nt2=2044; 

figure(1)
plot(Ubr)
title('measured bedload with MPM')
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 