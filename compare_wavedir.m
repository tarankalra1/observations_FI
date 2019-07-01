 clear all ; close all ; clc;  

 nt1=1900; 
 nt2=2044; 
 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
%
 netcdf_load(wh)
 Hs(:)=squeeze(wh_4061(1,1,:));
 Td(:)=squeeze(wp_peak(1,1,:));
 h(:)=squeeze(hght_18(1,1,:));
 Dwave_whp(:)=squeeze(wvdir(1,1,:)); 
% Dwave_whp=isnan(Dwave_whp);
% plot(squeeze(wvdir(1,1,nt1:nt2)))
%plot(Dwave_whp(nt1:nt2))
 %Dwave_whp(isnan(Dwave_whp))=[] ;
 vec=Dwave_whp;
 vec = vec(~isnan(vec))

 plot( vec(nt1:nt2) )
% % 
load('matfiles/skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','omega_br','Ubr','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
% 
% % for i=1:length(Dwave_whp)  
% %    if(~isnan(Dwave_whp(i)))
% %      if (Hs(i)>100);
% %        Hs(i)=0.0;
% %      Dwave_whp(i)=0.0; 
% %    end
% %     if (Td(i)>30);
% %         Td(i)=0.0;
% %     end
% %     [ubr_linear_wh(i),Tbav(i)]=ubspecfun( Hs(i),Td(i),h(i) );
% % end
% % 
% % hold on
% nt1=2000; nt2=2044; 
% plot(dn_wh(nt1:nt2),Dwave_whp(nt1:nt2),'linewidth',2);
%   hold on
%   plot(ang_rot(nt1:nt2)-180)
 %  plot(dn(nt1:nt2),ang_rot(nt1:nt2),'linewidth',2);
%  xlim([dn(nt1) dn(nt2)]);
%  datetick('x',2) % 'keeplimits')
% legend('Dwave-work horse','workhorse-pca')
% ylabel('Dwave (degrees)')
% print -dpng 'ubr_plot.png')
% 

