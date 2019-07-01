clear all ; close all ; clc; 

load('utot_art_waveform.mat','u_tot')

%plot(u_tot)
load('../matfiles/skewness_orbital_array.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
            'Hrmsu','Ubr','Ursell','dn','jtb_rec') 
       
% plot(ur_maj_rot_array(680,:))
% hold on
% plot(u_tot(29200:29200+40))