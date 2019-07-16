%% puv_proc_FI - Process ADV records from Fire Island ADV
% based on puv_proc_MVCO2011.m version of Feb 13, 2009
% last revised 5/1/2019
clear all ; close all ; clc ; 

kN1 = 30*0.018e-2 % roughness to be used on GM model
kN2 = 30*0.4e-2

% advbfn = '9885advb-cal.nc'; % burstfile name
% advsfn = '9885advs-cal.nc'; % statistics filename


advbfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advb-cal.nc'); % burstfile name
advsfn = fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9917advs-cal.nc'); % statistics filename

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

fs = ncreadatt(advbfn,'/','ADVDeploymentSetupSampleRate')
nsamp = ncreadatt(advbfn,'/','ADVDeploymentSetupSamplesPerBurst')
nominal_depth = ncreadatt(advbfn,'/','WATER_DEPTH') % nominal

%% overview plot
%figure(1);clf
%h1=plot(dn,u_1205,'linewidth',2)
%hold on
%h2=plot(dn,v_1206,'linewidth',2)
%xlim([dn(1) dn(end)]);
%datetick('x','keeplimits')
%ylabel('Velocity [cm/s]')
%ax1 = get(gca);
%ax1_pos = ax1.Position; % position of first axes
%ax2 = axes('Position',ax1_pos,...
%   'XAxisLocation','top',...
%   'YAxisLocation','right',...
%   'Color','none');
%bticks = [0:200:length(dn)]
%set(ax2,'Xticklabel',num2str(bticks'))
%set(ax2,'Xtick',bticks)
%set(ax2,'Ytick',[])
%set(ax2,'xlim',[bticks(1) bticks(end)])

% plot(dn,sqrt(u_1205.^2+v_1206.^2))
%figure(2); clf
%[sd1 az1 sd2 az2]=pcastats(u_1205,v_1206,25,1)

%% process bursts with no QA/QC
for n = 1:length(dn)
%  nt1=680; nt2=730; 
% count=1; 
% for n=nt1:nt2
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
      
      %u=detrend(u);
      %v=detrend(v); 
      % TODO - QA/QC, replace sketchy values here
      
      % TODO - Do we want to do any filtering here?
      
      % quick look at raw data
 
     figure(4); clf
     [sd1 az1 sd2 az2]=pcastats(u*100,v*100,50,1);
     % size(fs)
     % Get the statistics directly 
     UBS(n) = ubstatsr( u, v, fs );
     
     % get ubr, Hrms
     PUV(n) = puvq(p, (u), (v), depth(n), zp(n), zr(n), fs, 1050, 1030., 0.04, 1/6);
     
     if 1
     figure(5); clf
     plot(UBS(n).ur,UBS(n).vr,'.')
     title(sprintf('Rotated Burst Velocities burst = %d , %s',bn,dnsb))
     xlabel('ur, [m ^. s^{-1}]')
     ylabel('vr, [m ^. s^{-1}]')
     set(gca,'xlim',[-0.5 0.5]);
     set(gca,'ylim',[-0.5 0.5])
     pause(0.1)
     end
     
     kh = qkhfs( 2*pi/PUV(n).Tr, depth(n) );
     Tr(n)=PUV(n).Tr; 
     Ubr(n)=PUV(n).ubr ; 
     Hrmsu(n)=PUV(n).Hrmsu; 
     k(n) = kh./depth(n);
     Ur(n) = 0.75*0.5*PUV(n).Hrmsu*k(n)./(kh.^3); % RRvR Eqn. 6.
     %count=count+1
     %dnsb_rec(n)=dnsb ;
     jtb_rec(n)=jtb ;
     
     velu_skew(n)=UBS(n).ur_sk;
     velv_skew(n)=UBS(n).vr_sk; 
     save('crs_adv_9917.mat','dn','jtb_rec','depth','Hrmsu',..........
          'Tr','Ubr','Ur','velu_skew','velv_skew')
%       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Use the Nortek wds routines to check.
      % These are slightly modified from 
      % http://www.nortekusa.com/usa/knowledge-center/table-of-contents/waves
     lf=.035;      %Hz - low frequency cutoff
     maxfac=200;   %   - maximum value of factor scaling pressure to waves
     minspec=0.1;  %m^2/Hz - minimum spectral level for computing
                   %         direction and spreading
     Ndir=0;        %deg - direction offset (includes compass error and 
             %      misalignment of cable probe relative to case
   %           the offset for the Aquadopp Profiler is 0
%
     parms=[lf maxfac minspec Ndir];
     hp = nanmedian(zp);
     hv = -pdelz;
     nF = 1050;
     
   %  u=detrend(u); v=detrend(v); p=detrend(p);
     
  %   u=iwavesbp(u,0.04) ; %0.04 Hz 
     [Su,Sp,Dir,Spread,F,dF,DOF] = wds(u,v,p,1/fs,nF,hp,hv,parms);
     
     
     [Hs(n),peakF(n),peakDir(n),peakSpread(n)] = hs(Su,Sp,Dir,Spread,F,dF);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %save('nortek.mat','Hs', 
      
   end
end
%% The directions from PUV need to be flipped 180, I think.
azr = 180+[PUV(:).azr];

% and the negative directions from hs need to be converted to 0 - 360.
 peakDir(peakDir<0)=peakDir(peakDir<0)+360.;

%% make some summary plots
ok = find(~isnan(depth));
rp = ruessink_asymm( Ur );
Su = rp.Su;
Au = rp.Au;
r = rp.r;
sk = UBS.ur_sk;

save mat\puv_proc_FI
% 
% figure(5); clf
% subplot(411)
% h1=plot(dn(ok),[PUV(:).Hrmsu],'linewidth',2);
% hold on
% h2=plot(dn(ok),[PUV(:).Hrmsp],'.');
% h3=plot(dn,Hs,'-');
% xlim([dn(1) dn(end)]);
% datetick('x','keeplimits')
% ylabel('Hrmsu, Hrmsp [m]')
% legend([h1;h2;h3],'Hrmsu','Hrmsp','Hs wds')
% 
% 
% subplot(412)
% h1=plot(dn(ok),azr,'linewidth',2);
% hold on
% h2=plot(dn,peakDir,'.');
% h3=plot(dn(ok),180+[UBS(:).maj_az],'.');
% xlim([dn(1) dn(end)]);
% datetick('x','keeplimits')
% ylabel('Direction [\circT]')
% legend([h1;h2;h3],'Dp PUV + 180','Dp wds','PCA maj + 180')
% 
% subplot(413)
% plot(dn(ok),log10(Ur(ok)))
% xlim([dn(1) dn(end)]);
% ylabel('log10[ Ur ]')
% datetick('x','keeplimits')
% 
% subplot(414)
% plot(dn(ok),r(ok))
% hold on
% plot(dn(ok),Su(ok))
% % plot(dn(ok),[UBS(:).ur_sk])
% ylabel('r, Su')
% 
% xlim([dn(1) dn(end)]);
% datetick('x','keeplimits')
% print -dpng 'summary_plot.png') 