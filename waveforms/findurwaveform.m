function wf=findurwaveform(ur,fs,dnst)
%function wf=findurwaveform(ur,fs,dnst)
%script to find representative crest and trough wave form from nearbed
%velicity time series rotated into along wave direction. signal maybe
%bandpass filter to included only incident waves prior to rotation
%
%INPUT:
%ur = rotated, along wave, nearbed velocity burst 1xN (m/s)
%fs = sampling rate of burst data (Hz)
%dnst = burst start time matlab datenum

%
%OUTPUT:
%wf - sturct of individual waveforms in burst
%   .dn = wave form datenum start time
%   .t = wave form time (s)
%   .ub = wave nearbed velocity
%   .Tc = time crest (s)
%   .Tt = time trough (s)
%   .umax = amplitude crest (m/s)
%   .umin = amplitude trough (m/s)
%   .Ac = area crest (m/s * s)
%   .At = area trough (m/s *s)

% clear
% load('mat\puv_proc_FI_iwaves.mat')
% nb=149
% fs=8; %samp rate of signal (Hz)
% ur=[UBS(nb).ur]'; %along wave rotated, iwaves bandpass filtered near-bed veloicty (1xN)

t=[1:length(ur)]/fs; %burst time (s)

try
    %find zero crossing times
    tz=findzs(ur,t);


    %loop to find crest/trough info
    for ii=1:length(tz)-2
        %find half wave signal
        idx=find(t>=tz(ii) & t<=tz(ii+1));
        s(ii)={[0 ur(idx) 0]};
        ts(ii)={[tz(ii) t(idx) tz(ii+1)]};
        A(ii)=trapz(ts{ii},s{ii});
        S_max(ii)=max(s{ii});
        S_min(ii)=min(s{ii});
        T(ii)=tz(ii+1)-tz(ii);
    end

    cnt=0
    xlim=[0 15]
    if sign(A(1))>=0
        for i=1:2:length(ts)-1
        cnt=cnt+1
        wf(cnt).dn=dnst+[ts{i} ts{i+1}(2:end)]/(3600*24);
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)]
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).umax=S_max(i);
        wf(cnt).umin=S_min(i+1);
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).umax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).umin)-wf(cnt).Tc;
    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end

    else
        for i=2:2:length(ts)-2
        cnt=cnt+1
        wf(cnt).dn=dnst+[ts{i} ts{i+1}(2:end)]/(3600*24);
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)]
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).umax=S_max(i);
        wf(cnt).umin=S_min(i+1);
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).umax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).umin)-wf(cnt).Tc;

    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end
    end
catch
        wf.dn=dnst;
        wf.t=NaN;
        wf.ub=NaN;
        wf.T=NaN;
        wf.Tc=NaN;
        wf.Tt=NaN;
        wf.umax=NaN;
        wf.umin=NaN;
        wf.Ac=NaN;
        wf.At=NaN;
        wf.Tcu=NaN;
        wf.Ttu=NaN;
end

return