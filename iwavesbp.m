function sf=iwavesbp(s,fs,T_long,T_short)
%function sf=iwavesbp(s,fs)
%
%m-fcn to bandpass filter for incident waves if the 2-25s period range
%
%INPUT
% s = signal input (u,v,w, or press);
% fs = sampling frequency (Hz)
% T_long = longest wave period (s)
% T_short = shortest wave period (s)
%
%OUTPUT
% sf - bandpass filtered signal (same units as input)

N= 2; %filter order
Ny=fs/2; %Nyquist Freq (Hz)

Wn=[1/(T_long*Ny) 1/(T_short*Ny)]; %bandpass frequencies normalized by Nyquist

[b,a]=butter(N,Wn);

sf=filtfilt(b,a,s);
