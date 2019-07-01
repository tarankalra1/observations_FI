function sf=iwavesbp(s,fs,t_up,t_low)
%function sf=iwavesbp(s,fs)
%
%m-fcn to bandpass filter for incident waves if the 2-25s period range
%
%INPUT
% s = signal input (u,v,w, or press);
% fs = sampling frequency (Hz)
%
%OUTPUT
% sf - bandpass filtered signal (same units as input)

N= 2; %filter order
Ny=fs/2; %Nyquist Freq (Hz)

Wn=[1/(t_up*Ny) 1/(t_low*Ny)]; %bandpass frequencies normalized by Nyquist

[b,a]=butter(N,Wn);

sf=filtfilt(b,a,s);
