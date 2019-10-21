clear all ; close all ; clc ;
% Tarandeep S Kalra 
% Compare skewness with workhorse and directly from ADV 

nt1=1; nt2= 2044; 
% 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')

Su_skewness_old=0.0;     
for t=nt1:nt2
  Su_skewness_adv(t)=Su_skewness(t); %+Su_skewness_old; 
  Au_skewness_adv(t)=Au_skewness(t); %+Su_skewness_old; 
end 

% WORKHORSE DATA from linear wave theory
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)
 %Hs(:)=squeeze(wh_4061(1,1,:));    % extract significant wave  height 
  Td(:)=squeeze(wp_peak(1,1,:));    % extract peak wave period 
 %depth(:)=squeeze(hght_18(1,1,:)); % extract depth; 
% Dwave(:)=squeeze(wvdir(1,1,:)); 

% % Vspecdat 
 load('/media/taran/DATADRIVE2/Obs_data/matfiles/workhorse_emp_waveform_ubspecdat_vspec_inithght.mat',.....
     'Ur_emp','Hs','Tbr','h')
%       'umax_emp','umin_emp','Tc_emp','Tt_emp',........
%       'Tcu_emp','Ttu_emp','RR_emp','beta_emp','uhat_emp');
% %
% CALCULATE SKEWNESS FROM SURFACE WAVES using vspec
for i=nt1:nt2
  depth(i)=h(i) ; 
  omega=2.0*pi/Tbr(i);
  k=qkhfs(omega,depth(i))/depth(i);
  a_w=0.5*Hs(i);
  Ur(i)=0.75*a_w*k/((k*depth(i))^3.0);   

   
% Ursell number calcualtion from Ursell number and empirical relationships
% to get Su, Au..
  rp = taran_ruessink_empirical_skewness(Ur(i));
  Su_skewness_ruess(i)=rp.Su;
  Au_skewness_ruess(i)=rp.Au;
end   
 
%jt = time+time2/(3600*24*1000);
%dn_wh = j2dn(time,time2);
%dnsb_wh = datestr(datenum((dn_wh)));
  
% h=depth; % Depth from ADV  
% %ntime=end 
%for i=nt1:nt2
%  if (Hs(i)>100);
%    Hs(i)=0.0;
%  end
%  if (Td(i)>30); 
%    Td(i)=0.0;
%  end 

%end 

% figure(1)
% plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'k--');
% hold on 
% plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'r'); 
% legend('Direct','Parameterized')
% xlim([dn(nt1) dn(nt2)]);
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % ' 
% ylim([-0.3 0.3])
%  %print -dpng '-r300' '../pngfiles/Su_skewness.png'  
% % 
% figure(2)
%plot(dn(nt1:nt2),10*Su_skewness_adv(nt1:nt2));
%hold on
%ylabel('10*Su skewness adv')
%plot(dn(nt1:nt2),(Hrmsu(nt1:nt2)));
%plot(Td,Tbr)
%yyaxis right
%ylabel('Hrmsu')
%hold on 
%plot(dn(nt1:nt2),Hrmsu(nt1:nt2),'r--','linewidth',2); 
%legend('skewness adv-measured','Hrmsu')
%xlim([dn(nt1) dn(nt2)]);
%datetick('x',2) % 'keeplimits')
%print -dpng '../pngfiles/Su_skewness_Hrmsu.png' 


iremove=1 
  if(iremove==1)
% % 
% %
% figure(2) 
% plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'k--');
% hold on 
% plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'r'); 
% hold on 
% %hold on 
% plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2)*0.0,'linewidth',2); 
% legend('skewness adv-measured','skewness-ruessink','Location','Southwest')
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2,'keepticks','keeplimits');
% ylim([-0.3 0.3])


% REMOVE HRMSU <0.5 
for i=nt1:nt2
  if(Hs(i)<1.5)
    Su_skewness_adv(i)=0.0;
    Su_skewness_ruess(i)=0.0; 
   else 
    %Su_skewness_adv(i)=0.0;
    %Su_skewness_ruess(i)=0.0; 
    %Hrmsu(i)=0.0;
  end 
end
% 
% 
% % ARTIFICIALLY MAKE Sskewness to positive numbers
% Su_skewness_adv_abs=abs(Su_skewness_adv);

figure(3)
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2),'k');
hold on 
plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'r','LineWidth',2); 
hold on 
plot(dn(nt1:nt2),Su_skewness_adv(nt1:nt2)*0.0,'k--') ;%'linewidth',2); 
legend('Direct','Parameterized','Location','Northeast')
xlim([dn(nt1) dn(nt2)]);
datetick('x',2,'keepticks','keeplimits');
 
ylim([-0.22 0.22])

ylabel('S_{u}')
xlabel('Time (year/month/day)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig9_Su_skewness.png')
%print -dpng  '-r300' '../pngfiles/Su_skewness_Hrmsu_gt_0_7.png'  
end



% 
% xx=Su_skewness_adv_abs(nt1:nt2); 
% yy=Su_skewness_ruess(nt1:nt2)  ; 
% 
% figure(3)
% scatter(Su_skewness_adv_abs(nt1:nt2),Su_skewness_ruess(nt1:nt2),'r.','linewidth',2);
% R=corrcoef(xx,yy);
% RR=R(1,2); 
% title(['correlation coefficient is ',num2str(RR)]) %hold on 
% xlabel('Su skewness- ADV') ;
% ylabel('Su skewness-Ruessink calculated') 
% %scatter(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'linewidth',2);  
% %legend('skewness adv-measured','skewness-ruessink')
% %xlim([dn(nt1) dn(nt2)]);
% %datetick('x',2) % 'keeplimits')
% print -dpng '../pngfiles/Su_skewness_scatter.png'  
%  
