% written by Tarandeep S Kalra with help from Steve Suttles
% it calls the output of ADV based Hrmsu, ubr, Tr and compares
% with the workhorse data to get the work horse based wave orbital velocity and representative
% time period 
% compare pspec, vspec, dspec and ADV 

clear all ; close all ; clc; 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% convert work horse surface data of wave energy spectra to get ubr and
% Tbr..

netcdf_load(wh)
nt1=1; nt2=2044; 

Hs(:)=squeeze(wh_4061(1,1,:));
h_1d(:)=double(hght_18(1,1,:));
wdir(:)=squeeze(wvdir(1,1,:)); 

% Do we need this really ? 
initial_sensor_height = 2.11;

% Add the height to depth
h_1d(:)=h_1d(:)+initial_sensor_height ;


s_pspec(:,:)=double(pspec(1,1,:,:));  % pspec 
s_vspec(:,:)=double(vspec(1,1,:,:));  % vspec
s_dspec(:,:)=double(sspec(1,1,:,:));  %dspec

band_width=0.015625  ;
f=double(frequency(:,1));
df=band_width ;

isave=0 ; 

s_pspec=(s_pspec*0.001).^2 ; % converted to sq.m/Hz from sq.mm/Hz 
s_vspec=(s_vspec*0.001).^2 ; 
s_dspec=(s_dspec*0.001).^2 ;

for t=nt1:nt2; 
    
     %Hs_contour=Hs; 
     if (Hs(t)>100);
        Hs(t)=0.0;
     end
     Hs_contour(t)=Hs(t) ; 
      if (wdir(t)>1000); 
       wdir(t)=0.0; 
      end 
     wdir_contour(t)=wdir(t) ; 
     
     s_pspec(s_pspec<0)=0.0;
     s_pspec(s_pspec>1e12)=0.0; 
     
     s_vspec(s_vspec<0)=0.0; 
     s_vspec(s_vspec>1e12)=0.0;
     
     s_dspec(s_dspec<0)=0.0; 
     s_dspec(s_dspec>1e12)=0.0 ; 
     
     [uhat_pspec(t),Tbr_pspec(t)]=ubspecdat(squeeze(h_1d(t)),s_pspec(:,t)',f(:,1)',df);
     [uhat_vspec(t),Tbr_vspec(t)]=ubspecdat(squeeze(h_1d(t)),s_vspec(:,t)',f(:,1)',df); 
     [uhat_dspec(t),Tbr_dspec(t)]=ubspecdat(squeeze(h_1d(t)),s_dspec(:,t)',f(:,1)',df); 

end 
 


% ADV DATA 
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
Tbr_adv=Tr; 
%Hs_adv=sqrt(2.0)*Hrmsu; 
%Hs_adv=Hrmsu; 

%save('/media/taran/DATADRIVE2/Obs_data/matfiles/pspec_uhat_tr.mat','uhat_wh','Tr_wh')

isave1=0;
if(isave1==1)
figure(1)
plot(dn_wh(nt1:nt2),uhat_adv(nt1:nt2),'k--');
hold on
plot(dn_wh(nt1:nt2),uhat_pspec(nt1:nt2),'r');
hold on
plot(dn_wh(nt1:nt2),uhat_vspec(nt1:nt2),'b'); 
hold on
plot(dn_wh(nt1:nt2),uhat_dspec(nt1:nt2),'g') ; 
xlim([dn_wh(nt1) dn_wh(nt2)]);

datetick('x',2,'keepticks','keeplimits');
startdatenum=datenum('02-02-2014','mm-dd-yyyy');
enddatenum=datenum('05-04-2014','mm-dd-yyyy');
%set(gca,'XTick',[startdatenum:4:enddatenum]); 
legend('Direct','Pspec','Vspec','Sspec')
xlabel('Time (year/month/day)')
ylabel('Representative wave orbital velocity (m/s)')
%print -dpng '../pngfiles/repres/uhat_adv_wh_vspecdat.png'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig2_ubr.png')

figure(2)
scatter(uhat_adv,uhat_pspec,'r.')
R=corrcoef(uhat_adv,uhat_pspec,'rows','complete');
RR1=R(1,2);
hold on
scatter(uhat_adv,uhat_vspec,'b.')
R=corrcoef(uhat_adv,uhat_vspec,'rows','complete');
RR2=R(1,2);
hold on
scatter(uhat_adv,uhat_dspec,'g.')
R=corrcoef(uhat_adv,uhat_dspec,'rows','complete');
RR3=R(1,2);
hold on
plot(uhat_adv,uhat_adv,'k')
legend(sprintf('pspec, r^{2}=%g', RR1),sprintf('vspec, r^{2}=%g', RR2),..........
       sprintf('sspec, r^{2}=%g', RR3), 'Location', 'Northwest')
xlabel('Direct representative wave orbital velocity (m/s)')
ylabel('Wave spectra based representative wave orbital velocity (m/s)')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig2_ubr_scatter.png')

end 

isave2=0;
if(isave2==1)
figure(3)
plot(dn_wh(nt1:nt2),Tbr_adv(nt1:nt2),'k--');
hold on
plot(dn_wh(nt1:nt2),Tbr_pspec(nt1:nt2),'r');
hold on
plot(dn_wh(nt1:nt2),Tbr_vspec(nt1:nt2),'b'); 
hold on
plot(dn_wh(nt1:nt2),Tbr_dspec(nt1:nt2),'g') ; 
xlim([dn_wh(nt1) dn_wh(nt2)]);

datetick('x',2,'keepticks','keeplimits');
startdatenum=datenum('02-02-2014','mm-dd-yyyy');
enddatenum=datenum('05-04-2014','mm-dd-yyyy');
%set(gca,'XTick',[startdatenum:12:enddatenum]); 
legend('Direct','Pspec','Vspec','Sspec')
xlabel('Time (year/month/day)')
ylabel('Representative wave period (s)')
%print -dpng '../pngfiles/repres/uhat_adv_wh_vspecdat.png'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig2_Tbr.png')

figure(4)
scatter(Tbr_adv,Tbr_pspec,'r.')
R=corrcoef(Tbr_adv,Tbr_pspec,'rows','complete');
RR1=R(1,2);
hold on
scatter(Tbr_adv,Tbr_vspec,'b.')
R=corrcoef(Tbr_adv,Tbr_vspec,'rows','complete');
RR2=R(1,2);
hold on
scatter(Tbr_adv,Tbr_dspec,'g.')
R=corrcoef(Tbr_adv,Tbr_dspec,'rows','complete');
RR3=R(1,2);
hold on
plot(Tbr_adv,Tbr_adv,'k')
%text(8.0,12.0,['r^{2}_{pspec} = ' (num2str(RR1)) ])%,['r^{2}_{vspec} = ' (num2str(RR2)) ])

xlim([6 14])
ylim([6 14])
legend(sprintf('pspec, r^{2}=%g', RR1),sprintf('vspec, r^{2}=%g', RR2),..........
       sprintf('sspec, r^{2}=%g', RR3), 'Location', 'Northwest')   % r^{2}='num2str(RR1),'Vspec','Sspec','Location','Northwest')
xlabel('Direct representative wave period (s)')
ylabel('Wave spectra based representative wave period (s)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
print('-dpng','-r300','fig2_Tbr_scatter.png')

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6])
%print('-dpng','-r300','fig2_tbr.png')

end 

isave3=0;
if(isave3==1)
%   subplot(3,1,1)
%   pcolorjw(f,Dp,pspec)
%   subplot(3,1,2)
%   pcolorjw(f,Df,vspec)
%   subplot(3,1,3)
%   pcolorjw(f,D,dspec)
 
  figure(5)
  plot( (squeeze(pspec(1,1,:,155))*0.001).^2 )
  hold on 
  plot( (squeeze(vspec(1,1,:,155))*0.001).^2 )
  hold on;
  plot( (squeeze(sspec(1,1,:,155))*0.001).^2,'o' )
  ylim([0 30])
  legend('Pspec','Vspec','Sspec','Location','Northeast')
  xlabel('Frequency (Hz)')
  ylabel('Non-directional Wave Height Spectrum (m^{2}/Hz)')
  
  figure(6)
  plot( (squeeze(pspec(1,1,:,155)) ) ) 
  hold on 
  plot( (squeeze(vspec(1,1,:,155)) ) ) 
  hold on;
  plot( (squeeze(sspec(1,1,:,155)) ) ) 
  %ylim([0 30])
  legend('Pspec','Vspec','Dspec','Location','Northeast')
  xlabel('Frequency (Hz)')
  ylabel('Non-directional Wave Height Spectrum (m^{2}/Hz)')
    
end 

figure(7)
% ;
cdir=wdir_contour ;  
scatter(Tbr_adv,Tbr_vspec,25,cdir,'filled')
hold on 
plot(Tbr_adv,Tbr_adv,'k')
xlim([6 14])
ylim([6 14])
caxis([180 270])
colorbar
%legend('Pspec','Vspec','Sspec','Location','Northwest')
xlabel('Direct representative wave period (s)')
ylabel('Wave spectra based representative wave period (s)')
print('-dpng','-r300','fig2_Tbr_scatter.png')

%legend('num2str(RR1)','num2str(RR2']
% figure(2)
% scatter(uhat_adv,uhat_wh,'r.')
% R=corrcoef(uhat_adv,uhat_wh);
% RR=R(1,2);
% %title('Tr')
% title(['uhat, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
%legend('num2str(RR1)','num2str(RR2']
% figure(2)
% scatter(uhat_adv,uhat_wh,'r.')
% R=corrcoef(uhat_adv,uhat_wh);
% RR=R(1,2);
% %title('Tr')
% title(['uhat, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
%print -dpng '../pngfiles/repres/scatter_uhat_adv_wh_vspecdat.png'

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
% 
% figure(5)
% plot(dn_wh(nt1:nt2),Tr_adv(nt1:nt2),'k');
% hold on
% plot(dn_wh(nt1:nt2),Tr_wh(nt1:nt2),'r--');
% xlim([dn_wh(nt1) dn_wh(nt2)]);
% datetick('x',2) % 'keeplimits')% figure(2)
% legend('adv-Hs','workhorse-Hs')
% title('Tr')
% print -dpng '../pngfiles/repres/Tr_adv_wh_vspecdat.png'
%  
% figure(6)
% scatter(Tr_adv,Tr_wh,'r.')
%  R=corrcoef(Tr_adv,Tr_wh,'rows','complete');
%  RR=R(1,2);
% %title('Tr')
%  title(['Tr, correlation coeff is ',num2str(RR)])
% xlabel('ADV')
% ylabel('WH')
% print -dpng '../pngfiles/repres/scatter_Tr_adv_wh_vspecdat.png'
 