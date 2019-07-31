% written by Tarandeep S Kalra with help from Steve Suttles
% it calls the output of ADV based Hrmsu, ubr, Tr and compares
% with the workhorse data to get the work horse based wave orbital velocity and representative
% time period 

clear all ; close all ; clc; 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% convert work horse surface data of wave energy spectra to get ubr and
% Tbr..

netcdf_load(wh)
nt1=1; nt2=2044; 

Hs(:)=squeeze(wh_4061(1,1,:));
h_1d(:)=double(hght_18(1,1,:));
s(:,:)=double(vspec(1,1,:,:));  

band_width=0.015625  ;
f=double(frequency(:,1));
df=band_width ;

isave=0 ; 

s=(s*0.001).^2; % converted to sq.m/Hz from sq.mm/Hz 

for t=nt1:nt2; 
 % for ft=1:length(f)
  %if(s<-1e8)
     s(s<0)=0.0;
     s(s>1e12)=0.0;
  %end 
    [ubr(t),Tbr(t)]=ubspecdat(squeeze(h_1d(t)),s(:,t)',f(:,1)',df);
    Hs_wh(t)=Hs(t);   
    Hs_wh(Hs_wh>100)=0.0;
    %if (Hs_wh(i)>100);
    %    Hs(i)=0.0;
    %end
    
end 
uhat_wh=ubr;
Tr_wh=Tbr;
% To find correlation coefficient make Tr_wh=0 when Tr_wh=NaN;
%Tr_wh((Tr_wh)>35)=0.0 ; 
%Tr_wh(isnan(Tr_wh))=0;

%Hs_wh=sqrt(2.0)*Hrmsu; 

hold on 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Tr','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%            'Hrmsu','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
%dnsb_adv = datestr(datenum(gregorian(jtb_rec)))  ;
dn_adv=dn; 
dn_wh=dn_adv;
dnsb_adv = datestr(datenum(dn_adv))   ;

ok = find(~isnan(depth));
%ntime=30; 

uhat_adv=Ubr;
Tr_adv=Tr; 
Hs_adv=sqrt(2.0)*Hrmsu; 
%Hs_adv=Hrmsu; 

save('/media/taran/DATADRIVE2/Obs_data/matfiles/vspec_uhat_tr.mat','uhat_wh','Tr_wh')

isave=0
if(isave==1)
figure(1)
plot(dn_wh(nt1:nt2),uhat_adv(nt1:nt2),'k');
hold on
plot(dn_wh(nt1:nt2),uhat_wh(nt1:nt2),'r--');
xlim([dn_wh(nt1) dn_wh(nt2)]);
datetick('x',2) % 'keeplimits')% figure(2)
legend('adv-Uhat','workhorse-uhat')
print -dpng '../pngfiles/repres/uhat_adv_wh_vspecdat.png'

figure(2)
scatter(uhat_adv,uhat_wh,'r.')
R=corrcoef(uhat_adv,uhat_wh);
RR=R(1,2);
%title('Tr')
title(['uhat, correlation coeff is ',num2str(RR)])
xlabel('ADV')
ylabel('WH')
print -dpng '../pngfiles/repres/scatter_uhat_adv_wh_vspecdat.png'

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
% R=corrcoef(Hs_adv,Hs_wh);
% RR=R(1,2);
% title(['Hs, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
% print -dpng '../pngfiles/repres/scatter_hs.png'

figure(5)
plot(dn_wh(nt1:nt2),Tr_adv(nt1:nt2),'k');
hold on
plot(dn_wh(nt1:nt2),Tr_wh(nt1:nt2),'r--');
xlim([dn_wh(nt1) dn_wh(nt2)]);
datetick('x',2) % 'keeplimits')% figure(2)
legend('adv-Hs','workhorse-Hs')
title('Tr')
print -dpng '../pngfiles/repres/Tr_adv_wh_vspecdat.png'
 
figure(6)
scatter(Tr_adv,Tr_wh,'r.')
 R=corrcoef(Tr_adv,Tr_wh,'rows','complete');
 RR=R(1,2);
%title('Tr')
 title(['Tr, correlation coeff is ',num2str(RR)])
xlabel('ADV')
ylabel('WH')
print -dpng '../pngfiles/repres/scatter_Tr_adv_wh_vspecdat.png'
end
