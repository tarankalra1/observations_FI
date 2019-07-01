 clear all; close all ; clc;

url='9891whp-cal.nc' 
 
netcdf_load(url)
Hs(:)=squeeze(wh_4061(1,1,:));
Td(:)=squeeze(wp_peak(1,1,:)); 
h=11.5; 

%time_north=datenum(gregorian(double(time)+double(time2)/3600/24/1000));
%time_north=nctime2dn(time,time2)
% plot(time2)
% Filter the values of Hs, Td 
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
     [uhat(i),Tbav(i)]=ubspecfun(Hs(i),Td(i),h); 
 end


%plot(time_north)

%plot((Hs(1,1:end)))
%datetick('x',6)


%plot(Td, Tbav,'.')


%time_north=datenum(gregorian(double(time)+double(time2)/3600/24/1000));

%plot(time_north(1:3000),1:3000)
%datetick('x',6)
%for 
% TODO - Should this be calculated from Hs?
%uhat
%Tbav
%dn1= nctime2dn(double(time),double(time2));
% julian2datenum(time_north)
% % 
% igd=squeeze(Td)<30;
% plot(dn1(igd),squeeze(Td(igd)),'.')
% igd=squeeze(wvdir)<1000; 
% plot(dn1(igd),squeeze(wvdir(1,1,igd)),'.')
% datetick('x')
