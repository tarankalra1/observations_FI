clear all ; close all ; clc ;
% Get ubr from ADV and compare with the ubr from workhorse
%% LOAD ADV DATA
%load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/crs_adv_9917.mat') 
     %save('crs_adv_9885.mat','dn','jtb_rec','depth','Hrmsu','Tr','Ubr')
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
%dnsb_adv = datestr(datenum(gregorian(jtb_rec)))  ;
dn_adv=dn; 
dnsb_adv = datestr(datenum(dn_adv))   ;

ok = find(~isnan(depth));
%ntime=30; 
nt1=1; nt2=2044; 

uhat_adv=Ubr;
Tr_adv=2*pi./omega_br; 
Hs_adv=sqrt(2.0)*Hrmsu; 
%Hrms; 
 

%% WORKHORSE DATA to get uwave_rms 
wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)
 Hs(:)=squeeze(wh_4061(1,1,:));
 Td(:)=squeeze(wp_peak(1,1,:));
 h(:)=squeeze(hght_18(1,1,:)); 
 
jt = time+time2/(3600*24*1000);
dn_wh = j2dn(time,time2);
dnsb_wh= datestr(datenum((dn_wh)));

for i=nt1:nt2
%   if(~isnan((Hs(i))))
%        Hs(i)=0;
%        Td(i)=0;
%        h(i)=0;
%   else
%        Hs(i)=Hs(i); 
%        Td(i)=Td(i);
%        h(i)=h(i);
%   end
    if (Hs(i)>100);
       Hs(i)=0.0;
    end
    if (Td(i)>30); 
        Td(i)=0.0;
    end 
     Hs_wh(i)=Hs(i); 
     [uhat_wh(i),Tr_wh(i)]=ubspecfun( Hs(i),Td(i),h(i) ); 
end
%  
% figure(1)
% plot(dn_wh(nt1:nt2),uhat_adv(nt1:nt2),'k');
% hold on 
% plot(dn_wh(nt1:nt2),uhat_wh(nt1:nt2),'r--');
% xlim([dn_wh(nt1) dn_wh(nt2)]);
% datetick('x',2) % 'keeplimits')% figure(2)
% legend('adv-Uhat','workhorse-uhat')
% print -dpng '../pngfiles/repres/uhat_adv_wh.png'
% 
% 
% figure(2)
% scatter(uhat_adv,uhat_wh,'r.')
% R=corrcoef(uhat_adv,uhat_wh)
% RR=R(1,2)
% title(['uhat, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
% print -dpng '../pngfiles/repres/scatter_uhat_adv_wh.png'
% 
% 
% figure(3)
% plot(dn_wh(nt1:nt2),Hs_adv(nt1:nt2),'k');
% hold on
% plot(dn_wh(nt1:nt2),Hs_wh(nt1:nt2),'r--');
% xlim([dn_wh(nt1) dn_wh(nt2)]);
% datetick('x',2) % 'keeplimits')% figure(2)
% legend('adv-Hs','workhorse-Hs')
% print -dpng '../pngfiles/repres/hs_adv_wh.png'
% 
% figure(4)
% scatter(Hs_adv,Hs_wh,'r.')
% R=corrcoef(Hs_adv,Hs_wh)
% RR=R(1,2)
% title(['Hs, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
% print -dpng '../pngfiles/repres/scatter_hs_adv_wh.png'

figure(5)
plot(dn_wh(nt1:nt2),Tr_adv(nt1:nt2),'k');
hold on 
plot(dn_wh(nt1:nt2),Td(nt1:nt2),'r--');
xlim([dn_wh(nt1) dn_wh(nt2)]);
datetick('x',2) % 'keeplimits')% figure(2)
legend('adv-Hs','workhorse-Hs')
title('Tr')
%print -dpng '../pngfiles/repres/Tr_adv_wh.png'
% 
% figure(6)
% scatter(Tr_adv,Tr_wh,'r.')
% R=corrcoef(Tr_adv,Tr_wh)
% RR=R(1,2)
% title('Tr')
% %title(['Tr, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
% print -dpng '../pngfiles/repres/scatter_Tr_adv_wh.png'
