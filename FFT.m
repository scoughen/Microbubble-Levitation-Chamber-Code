%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script takes a signal times series and calculates the fast fourier
%  transform.
%
%  The parameters that need to be set are:
%    Scan parameter:
%      file = file name of the signal time series
%      L = length of signal (the number of data pts in the signal)
%      T = the sampling period
% 
%  S. Coughenour - Nov. 17, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

file = "DrivingSignal_DriveOnly_Data.csv";
S = readmatrix(file);
S = S(:,1)';
S2 = [zeros(1,length(S)) , S];

T = 2e-9;             % Sampling period       
Fs = 1/T;             % Sampling frequency                    
L = 350000;           % Length of signal
t = (0:L-1)*T;        % Time vector

% X = S + 2*randn(size(t));
% 
% figure
% plot(1000*t(1:50),X(1:50))
% title("Signal Corrupted with Zero-Mean Random Noise")
% xlabel("t (milliseconds)")
% ylabel("X(t)")
% 
% Y = fft(X);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
f = Fs*(0:(L/2))/L;
% figure
% plot(f,P1) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")


Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
stem(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

figure
stem(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")
xlim([0 1e5])







