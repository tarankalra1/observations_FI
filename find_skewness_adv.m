%% puv_proc_FI - Process ADV records from Fire Island ADV
% based on puv_proc_MVCO2011.m version of Feb 13, 2009
% last revised 5/1/2019
clear all ; close all ; clc ; 

kN1 = 30*0.018e-2 ; % roughness to be used on GM model
kN2 = 30*0.4e-2   ; 

% advbfn = '9885advb-cal.nc'; % burstfile name
% advsfn = '9885advs-cal.nc'; % statistics filename

% advbfn = fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9917advb-cal.nc'); % burstfile name
% advsfn = fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9917advs-cal.nc'); % statistics filename

advbfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advb-cal.nc'); % burstfile name
advsfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advs-cal.nc'); % statistics filename

%advbfn = fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9885advb-cal.nc');
%advsfn = fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9885advs-cal.nc'); 

ncload(advsfn); % load the statistics file

zr = 0.4; % placeholder...need to check measurement elevation here
bps.instrname = 'ADV 9917';
%bps.instrname= 'ADV 9885'; 

% you can see all of the variables in the burst file with:
% ncdisp(advbfn);
%time=ncread(advbfn,'time');
%time2=ncread(advbfn,'time2') ; 

jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

gb_first = 1; % assume first burst is good
gb_last = length(dn); % there are some bad bursts
gb = (gb_first:gb_last)';

dn = dn(gb);
brange = brange(gb); % Height of ADV transducer above boundary
zoff = ncreadatt(advbfn,'/','ADVProbeSamplingVolumeOffset')/100.
zr_median = nanmedian(brange)-zoff % this is the median ADV sample elevation
zr = brange-zoff; % time series of ADV sample locations

% zr has some pretty low values, indicating that measurements were being
% made 2 cmab. We need to check this and think about it.

z_init = ncreadatt(advbfn,'u_1205','initial_sensor_height')

ap = 10.13 % std atmos. pressure (or a time series from nearby) [dBar]
p_z = ncreadatt(advbfn,'P_4022','initial_sensor_height');
pdelz = p_z-z_init; % elevation diff. between velocity and pressure
zp = zr+pdelz; % elevation of pressure measurements [m] accounting for variable brange
depth = zr+pdelz+(P_4023(gb)*0.01-ap); % time series of depth [decibars ~= meters]
depth(depth>1e30)=NaN; %convert fill_values to NaNs
fs = ncreadatt(advbfn,'/','ADVDeploymentSetupSampleRate')
nsamp = ncreadatt(advbfn,'/','ADVDeploymentSetupSamplesPerBurst')
nominal_depth = ncreadatt(advbfn,'/','WATER_DEPTH') % nominal
 
%% process bursts with no QA/QC
%for = 1:length(dn)
nt1=140; nt2=180; 
%count=1; 
%nt2=nt1 ; 
%nt1=1;  nt2=2087; 
isave=1; 

% set the time period in seconds for excluding infragravity wave band
t_up=20; t_low=4 ;

for n=nt1:nt2  
   if(~isnan(depth(n)))
      bn = ncread(advbfn,'burst',n,1);      % this burst number from beginning...might just want to go from 1 to nb
      jtb = double(ncread(advbfn,'time',[1 n],[1 1]))+......
            double(ncread(advbfn,'time2',[1 n],[1 1])/(3600*24*1000)); 
      dnsb = datestr(datenum(gregorian(jtb))); 
      fprintf(1,'Burst %d at %s\n',bn,dnsb);
      u = ncread(advbfn,'u_1205',[1 n],[Inf 1])/100;
      v = ncread(advbfn,'v_1206',[1 n],[Inf 1])/100;
      w = ncread(advbfn,'w_1204',[1 n],[Inf 1])/100;
      p = ncread(advbfn,'P_4022',[1 n],[Inf 1])/100;
      a1 = ncread(advbfn,'AGC1_1221',[1 n],[Inf 1]);
      a2 = ncread(advbfn,'AGC2_1222',[1 n],[Inf 1]);
      a3 = ncread(advbfn,'AGC3_1223',[1 n],[Inf 1]);
      c1 = ncread(advbfn,'cor1_1285',[1 n],[Inf 1]);
      c2 = ncread(advbfn,'cor2_1286',[1 n],[Inf 1]);
      c3 = ncread(advbfn,'cor3_1287',[1 n],[Inf 1]);
      
      % DETREND U, V
      u_detrend=detrend(u); 
      v_detrend=detrend(v); 
      
      % BAND PASS FILTER ; 
      u_band=iwavesbp(u_detrend, fs, t_up, t_low); %0.04 Hz in seconds 
      v_band=iwavesbp(v_detrend, fs, t_up, t_low); 
      
      % MEDIAN FILTER
      u_med=medfilt(u_band, 9) ; %smooth
      v_med=medfilt(v_band, 9) ;
      
      u_send=u_band;
      v_send=v_band; 
        
%       % PCA STATS ;
        figure(1); clf
      [sd1 az1 sd2 az2]=pcastats(u_send*100,v_send*100,50,1);

      UBS(n) = ubstatsr( u_send, v_send, fs );
      ur_maj_rot=UBS(n).ur; % major rotated
      vr_min_rot=UBS(n).vr; % minor 
      ang_rot(n)=UBS(n).maj_az; 
      % 
%      FIND SKEWNESS 
      Su_skewness(n)=mean(ur_maj_rot.^3)/(std(ur_maj_rot)).^3;  % Eqn. 5 in the skewness calculation 
       
      hilbert_asym=imag(hilbert(ur_maj_rot)) ; 
      Au_skewness(n)=mean(hilbert_asym.^3)/(std(ur_maj_rot)).^3 ; 
      
      ur_bar(n)=mean(ur_maj_rot); 
      ur_cube(n)=mean(ur_maj_rot.^3); 
     % get ubr, Hrms
       PUV(n) = puvq(p, (u), (v), depth(n), zp(n), zr(n), fs, 1050, 1030., 0.04, 1/6);
      
       kh = qkhfs( 2*pi/PUV(n).Tr, depth(n) );
       Tr(n)=PUV(n).Tr;
       omega_br(n)=PUV(n).omegar ;  
       Ubr(n)=PUV(n).ubr ; 
       Hrmsu(n)=PUV(n).Hrmsu; 
     
       k(n) = kh./depth(n);
       Ursell(n) = 0.75*0.5*PUV(n).Hrmsu*k(n)./(kh.^3); % RRvR Eqn. 6.
     
       jtb_rec(n)=jtb ;
              
       if(isave==1)
         save('skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Tr','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
       end
       if(isave==2)
        
        ur_maj_rot_array(:,n)=ur_maj_rot; 
        vr_min_rot_array(:,n)=vr_min_rot; 
        save('skewness_orbital_array.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
           'Hrmsu','Ubr','Ursell','dn','jtb_rec') 
       end 
       
        elseif isnan(depth(n)) & exist('PUV')
            jtb = double(ncread(advbfn,'time',[1 n],[1 1]))+......
            double(ncread(advbfn,'time2',[1 n],[1 1])/(3600*24*1000));
    
        %set struct variblles to NaNs
            puvfnames=fieldnames(PUV);
                for i=1:length(puvfnames)
                PUV(n).(puvfnames{i})=nan(size(PUV(n-1).(puvfnames{i})));
                end 
           
            ubsfnames=fieldnames(UBS);
                for i=1:length(ubsfnames)
                UBS(n).(ubsfnames{i})=nan(size(UBS(n-1).(ubsfnames{i})));
                end
    
                %other index vars to NaN
                 Tr(n)=NaN; 
                 Ubr(n)=NaN; 
                 Hrmsu(n)=NaN; 
                 k(n) = NaN;
                 Ur(n) = NaN; 

                velu_skew(n)=NaN;
                velv_skew(n)=NaN;

                ang_rot(n)=NaN;
                 ur_maj_rot=UBS(n).ur; % major rotated
                 vr_min_rot=UBS(n).vr; % minor 
      
                Su_skewness(n)=NaN;    
                Au_skewness(n)=NaN; 
      
               ur_bar(n)=mean(ur_maj_rot); 
               ur_cube(n)=mean(ur_maj_rot.^3);  
                omega_br(n)=NaN ;  
                Ursell(n) = NaN;
     
                jtb_rec(n)=jtb ;
              
           if(isave==1)
             save('skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
                'Hrmsu','Tr','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
           end
           if(isave==2)

            ur_maj_rot_array(:,n)=ur_maj_rot; 
            vr_min_rot_array(:,n)=vr_min_rot; 
            save('skewness_orbital_array.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
               'Hrmsu','Ubr','Ursell','dn','jtb_rec') 
           end 
   end
end
%(nt1:nt2),Su_skewness_adv(nt1:nt2),'linewidth',2);
% hold on 
% plot(dn(nt1:nt2),Su_skewness_ruess(nt1:nt2),'linewidth',2); 
% 
% legend('skewness adv-measured','skewness-ruessink')
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % '
 %print -dpng 'Su_skewnes% figure(2)
% plot(UBS(n).ur*100, UBS(n).vr*100,'.')
%  %plot(ur,vr,'.')
% xlim([-50.0 50.0]) 
% ylim([-50.0 50.0])

figure(1)
plot(ur_cube)
%hold on 
%plot(u_med)










%legend('u','median')
% plot(dn
% %th = pi/4:pi/4:2*pi;
% th = linspace(0,2*pi,20);
% 
%figure(1)
% pax = gca;
% pax.ThetaDir = 'clockwise';
% 
% polarscatter(th,u,v)

% hold on
% plot(detrend(u))
% hold on
%u_band=iwavesbp(detrend(u), 25); %0.04 Hz in seconds 
%plot(u_band)
