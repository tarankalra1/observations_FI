
clear all ; close all ; clc ;
%  find gradient of cumbedload 

load('/media/taran/DATADRIVE2/Obs_data/matfiles/vandera_bedld_unet.mat',.....
    'bedldx_wh_unet','bedldy_wh_unet','R','Beta','Ur')
bedld_unet_time=bedldx_wh_unet; 
cumbedld_bedld_unet_time=cumtrapz(bedld_unet_time)*3600; 


grad_cumbed=bedld_unet_time(2:end)-bedld_unet_time(1:end-1); 

% get the current velocity magnitude
load('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_Steve/unet_9917adv_phi.mat',.......
                                               'phid','wdir','cdir','unet')


figure(1)
subplot(2,1,1)
plot(grad_cumbed(1220:1231),'ro--')
ylabel('chanage bedload')

subplot(2,1,2)
plot(cosd(phid(1220:1231)),'ro--')
ylabel('cosd phi')

