clear all ; close all ; clc ;
nt1=1 ; nt2= 2044; 
% all filters with different upper and lower frequency 
load('skewness_steve_4s_20sec.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_allfilter_4s_20s=Su_skewness; 

% all filters with different upper and lower frequency 
load('skewness_steve_4s_25sec.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_allfilter_4s_25s=Su_skewness; 

% without median filter (band pass + detrending)
load('skewness_steve_detrend_band_4_20.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_nomedian_4s_20s=Su_skewness; 
Ubr_nomedian=Ubr; 

% without median and band pass filter (detrending only)
load('skewness_steve_detrend_only.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_detrend=Su_skewness; 

% without median and band pass and detrend filter 
load('skewness_steve_nofilter.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_nofilter=Su_skewness;

% all median and band pass and detrend filter 
load('skewness_steve_allfilter_9.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
Su_skewness_allfilter_9=Su_skewness;
Ubr_allfilter_9=Ubr; 
% 
% % all median and band pass and detrend filter 
% load('skewness_steve_allfilter_med_9.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%             'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube','ang_rot','Au_skewness');
% Su_skewness_allfilter_9=Su_skewness;














figure(1)
plot(dn(nt1:nt2),Su_skewness_detrend(nt1:nt2),'ro-','linewidth',2);
hold on 
%plot(dn(nt1:nt2),Su_skewness_allfilter_4s_25s(nt1:nt2),'linewidth',2); 
%hold on 
plot(dn(nt1:nt2),Su_skewness_nomedian_4s_20s(nt1:nt2),'b--','linewidth',2); 
%hold on 
%plot(dn(nt1:nt2),Su_skewness_detrend(nt1:nt2),'linewidth',2); 
%hold on
%plot(dn(nt1:nt2),Su_skewness_nofilter(nt1:nt2),'linewidth',2); 
%legend('all filter 4-20s','all filter 4-25s','no median 4-20s', 'detrend only','no filter')
legend('all filter 4-20s','no median 4-20s')

xlim([dn(nt1) dn(nt2)]);
datetick('x',2) % '
% print -dpng 'skewness_with_withoutmedianfilter.png')   


% 
% figure(2)
% plot(dn(nt1:nt2),Ubr_allfilter_9(nt1:nt2),'ro-','linewidth',2);
% hold on 
% %plot(dn(nt1:nt2),Su_skewness_allfilter_4s_25s(nt1:nt2),'linewidth',2); 
% %hold on 
% plot(dn(nt1:nt2),Ubr_nomedian(nt1:nt2),'b--','linewidth',2); 
% %hold on 
% %plot(dn(nt1:nt2),Su_skewness_detrend(nt1:nt2),'linewidth',2); 
% %hold on
% %plot(dn(nt1:nt2),Su_skewness_nofilter(nt1:nt2),'linewidth',2); 
% %legend('all filter 4-20s','all filter 4-25s','no median 4-20s', 'detrend only','no filter')
% legend('all filter 4-20s','no median 4-20s')
% 
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % '
% print -dpng 'Ubr_with_withoutmedianfilter.png')   
%
