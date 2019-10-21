clear all ; close all ; clc ;
% Get ubr from ADV and compare with the ubr from workhorse
%% LOAD ADV DATA
%load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/crs_adv_9917.mat') 
     %save('crs_adv_9885.mat','dn','jtb_rec','depth','Hrmsu','Tr','Ubr')
% load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%             'Hrmsu','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
% %dnsb_adv = datestr(datenum(gregorian(jtb_rec)))  ;
% dn_adv=dn; 
% dnsb_adv = datestr(datenum(dn_adv))   ;
% 
% ok = find(~isnan(depth));
% %ntime=30; 
% nt1=1; nt2=2044; 
% 
% figure(1)
% plot(dn(nt1:nt2),Ubr(nt1:nt2),'linewidth',2); 
% xlim([dn(nt1) dn(nt2)]);
% datetick('x',2) % 'keeplimits')
%  
% %% Workhorse data that is using wave spectra to get uwave_rms
% load('/media/taran/DATADRIVE2/Obs_data/FI_processing_Tsk/matfiles/ubr_from_spectra.mat','ubr','Tbr')
% ubr_spectra_wh=ubr; 

%% WORKHORSE DATA to get uwave_rms 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)
 Hs(:)=squeeze(wh_4061(1,1,:));
 Td(:)=squeeze(wp_peak(1,1,:));
 h(:)=squeeze(hght_18(1,1,:)); 
 wdir(:)=squeeze(wvdir(1,1,:)); 
 
 nt1=1; nt2=2044; 

jt = time+time2/(3600*24*1000);
dn_wh = j2dn(time,time2);
dnsb_wh = datestr(datenum((dn_wh)));
% h=depth ; 
% h=isnan(depth); 
% h=depth; % Depth from ADV  
% %ntime=end 
count_norm=0; count=0 ;
count_per=0;  count_per_ex=0;
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
     if (wdir(i)>1000); 
         wdir(i)=0.0; 
     end 
     count_norm=count_norm+1; % count all ;
     if(Hs(i)>1.5);           % count waves > 1.5
        count=count+1;
     end 
     count_per=count_per+1; 
     if(Td(i)>10.0)
      count_per_ex=count_per_ex+1; 
     end 
         
     
    % [ubr_linear_wh(i),Tbav(i)]=ubspecfun( Hs(i),Td(i),h(i) ); 
 end
ifig=0; 
% figure(1)
 %set(gca,'fontWeight','bold','FontSize',12)
if(ifig==1)
figure(1)
set(gca, 'FontName', 'Times New Roman')
% subplot(3,1,1)
%set(gca,'FontSize', 24)
subplot(3,1,1)
plot(dn_wh(nt1:nt2),Hs(nt1:nt2),'r')
xlim([dn_wh(nt1) dn_wh(nt2)])
datetick('x',2,'keepticks','keeplimits')
ylabel('H_{s}(m)')
%xlabel('Time (year/month/day)')

subplot(3,1,2)
plot(dn_wh(nt1:nt2),Td(nt1:nt2),'r')
xlim([dn_wh(nt1) dn_wh(nt2)])
datetick('x',2,'keepticks','keeplimits')
ylabel('Peak wave period (s)')
%xlabel('Time (year/month/day)')

subplot(3,1,3)
plot(dn_wh(nt1:nt2),wdir(nt1:nt2),'r')
xlim([dn_wh(nt1) dn_wh(nt2)])
datetick('x',2,'keepticks','keeplimits');
startdatenum=datenum('02-02-2014','mm-dd-yyyy');
enddatenum=datenum('05-04-2014','mm-dd-yyyy');
%set(gca,'XTick',[startdatenum:12:enddatenum]);
% datetick('x',2,'keeplimits')
ylabel('Direction ({\circ})')
xlabel('Time (year/month/day)')

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 7])
print('-dpng','-r300','fig1.png')
end 

ifig=2 ;  
if(ifig==2)
figure(3)
WindRose(wdir+180,Td,'vwinds',[0 3 5 8 12],'TitleString',{''},'LabLegend','Peak wave period (s)','LegendVariable','T_{p}',.....
    'labels',{'North(0°)','South (180°)','East (90°)','West (270°)'})
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
 print('-dpng','fig1_dwave.png','-painters')
end 
ifig=3 ; 
if(ifig==3)
wind_rose(wdir,Td) %'vwinds',[0 3 5 8 12],'TitleString',{''},'LabLegend','Peak wave period (s)','LegendVariable','T_{p}',.....
end 
%    'labels',{'North(0°)','South (180°)','East (90°)','West (270°)'})
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
% print('-dpng','fig1_dwave.png','-painters')
 % 
% figure(4)
% WindRose(wdir,Hs,'vwinds',[0 0.5 1.0 1.5 3.0],'LabLegend','Significant wave height (s)','LegendVariable','H_{s}',.....
%     'labels',{'North(0°)','East (90°)','South (180°)','West (270°)'})

 % print('-dpng', 'WindRose.png','-painters');                                                         
%  hold on 
%  plot(dn_wh(nt1:nt2),ubr_linear_wh(nt1:nt2),'linewidth',2);
%  hold on 
%  plot(dn_wh(nt1:nt2),ubr_spectra_wh(nt1:nt2),'linewidth',2); 
%  xlim([dn_wh(nt1) dn_wh(nt2)]);
%   datetick('x',2) % 'keeplimits')
%   legend('adv-Ubr from burst data','workhorse-linear wave theory','workhorse-spectra')
%   ylabel('Ubr')
 %print -dpng 'ubr_plot.png'   
 