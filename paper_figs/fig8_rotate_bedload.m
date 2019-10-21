clear all ; close all ;clc;
%----Written by : Taran Kalra and Steve Suttles--------
%Use the vandera bedload that is aligned with the wave axis and perpendicular
%to it
%Use the shoreline angle and rotate the bedload along the shoreline and perpendicular
%to it to get cross-shore and along-shore bedload a
% find the angle between shoreline and dominant wave direction
% Rotate bedload to shoreline parallel and shoreline normal directions 
% \ /
%---------------------------------------

netwc=1
if(netwc==1) ; 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet.mat',.....
    'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur')
bdx_wdir_unet=bedldx_wh_unet; 
bdy_wdir_unet=bedldy_wh_unet;  
end 

wavesonly=1; 
%
if(wavesonly==1)
  load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_workhorse_vspec_inithght.mat',......
                                                           'bedldx_wh_vspec','R','Beta','Ur')
bdx_wdir=bedldx_wh_vspec;
bdy_wdir=bedldx_wh_vspec.*0.0; % cumtrapz(bedldx_empirical_vand_all);

end 

% Got FI shoreline 
shore_FI=importdata('/media/taran/DATADRIVE2/Obs_data/matfiles/shoreline_FI_strip.txt');
lon_shore_FI=shore_FI(:,1); lat_shore_FI=shore_FI(:,2);
shoreline_angle=shore_FI(:,3);
% Find mean of shoreline angle 
FI_angle=mean(shoreline_angle); 

% load the wave direction
 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/mat_9917adv_phi.mat','wdir')

nt1=1; nt2= 2044; 
load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_orbital_array.mat','dn') 
 
[uw,vw]=xycoord(1,wdir); %find wave velocities of unit vector in u,v coordinates
 
% Unit vector along x and y direction for shoreline 
[vs,us]=xycoord(1,FI_angle); % FI_angle is not polar but this function assumes polar, so used vs,us
%  
 
for i=1: length(wdir)
  if(~isnan(wdir(i)))  
    phir(i)=[ atan2(vw(i),uw(i)) - atan2(vs,us)  ]*180/pi         ;  
    bdx_sh(i)=cosd(phir(i))*bdx_wdir(i)-sind(phir(i))*bdy_wdir(i) ;
    bdy_sh(i)=sind(phir(i))*bdx_wdir(i)+cosd(phir(i))*bdy_wdir(i) ;
    
    bdx_wdir_unet_sh(i)=cosd(phir(i))*bedldx_wh_unet(i)-sind(phir(i))*bedldy_wh_unet(i); 
    bdy_wdir_unet_sh(i)=sind(phir(i))*bedldx_wh_unet(i)+cosd(phir(i))*bedldy_wh_unet(i); 
  end 
end 
 
dt=3600;
%  
figure(1) 

plot(dn(nt1:nt2),cumtrapz(bdx_sh)*dt,'k')
hold on
plot(dn(nt1:nt2),cumtrapz(bdx_wdir_unet_sh)*dt,'r')
xlim([dn(nt1) dn(nt2)]); 
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('Along-shore cummulative bedload, \int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
legend('Wave driven', 'Wave and current driven'); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print('-dpng','-r300','fig8_bdldx_unet_rotate.png')
% 
figure(2)
%subplot(2,1,2)

plot(dn(nt1:nt2),cumtrapz(bdy_sh)*dt,'k')
hold on 
plot(dn(nt1:nt2),cumtrapz(bdy_wdir_unet_sh)*dt,'r')
xlim([dn(nt1) dn(nt2)]); 
datetick('x',2,'keepticks','keeplimits');
xlabel('Time (year/month/day)')
ylabel('Cross-shore cummulative bedload, \int_{0}^{t} q_{b}\it{dt}\rm (m^{2})')
legend('Wave driven', 'Wave and current driven'); 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print('-dpng','-r300','fig8_bdldy_unet_rotate.png')