clear all ; close all ; clc; 

% CHECK WAVE ORBITAL VELOCITY 
 
 load('skewness_orbital_array.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec')  
nt1=680 ;nt2=730 ; 

% Positive skewness can be checked by trapz
z=trapz(ur_maj_rot_array(:,694)) 

% Negative skewness can be checked by trapz
z=trapz(ur_maj_rot_array(:,680)) 

% max at 694 
% % Plot positive skewness 
 figure(1)
plot(ur_maj_rot_array(:,694).*0.0)
hold on 
plot(ur_maj_rot_array(:,694))
title(' positive skewness')
ylabel('Time series of wave orbital velocity')
print -dpng 'positive_skewness.png'  

% trapz(

 figure(2)
% Plot negative skewness 
plot(ur_maj_rot_array(:,680).*0.0)
hold on 
plot(ur_maj_rot_array(:,680))
title('negative skewness')
ylabel('Time series of wave orbital velocity')
print -dpng 'negative_skewness.png'  

%xlim([dn(nt1) dn(nt2)]);
%datetick('x',2) % 'keeplimits')

%figure(3)
%plot(Su_skewness(680:730))