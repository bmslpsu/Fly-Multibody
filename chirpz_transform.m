function [Z,Fz,Mag,Phase] = chirpz_transform(t,x,f1,f2,n)
%% chirpz_transform: Computes DFT using Chirp Z-Transform
%
%   INPUT:
%       t   :  	time vector [s]
%       x   : 	data
%       f1  :   lower frequency bound
%       f2  :   upper frequency bound
%       n   :   # of frequency bins
%
%   OUTPUT:
%       Z   :   Chirp Z-Transform
%       Fz  :  	frequency vector [Hz]
%       Mag : 	magnitude vector
%       Phs : 	phase vector [rad]
%

% Sampling frequency [Hz]
Fs = 1 / mean(diff(t));

% Make frequency vector
fnz = (0:n-1)'/n;
Fz = (f2-f1)*fnz + f1;

% Compute Chirp-Z Transform
w = exp(-1j*2*pi*(f2-f1)/(n*Fs));
a = exp(1j*2*pi*f1/Fs);
Z = czt(x,n,w,a) ./ (n * Fs/100);

% Get magniyude & phase
Mag = 2*abs(Z);
Phase = angle(Z);

%% TESTING
% clear;close all ; clc
% 
% Fs = 1000;
% Fn = Fs/2;
% t = (0:(1/Fs):10)';
% h = 5*cos(2*pi*5*t) + 4*cos(2*pi*4*t) + + 3*cos(2*pi*3*t + pi);
% 
% L = length(h);
% n = L;
% y = fft(h) ./ n;
% fny = (0:n-1)'/n;
% % fy = Fs*fny;
% 
% fy = (linspace(0, 1, fix(L/2)+1)*Fn)';  % frequency vector [Hz]
% Iv = 1:length(fy);                  	% index vector
% y = y(Iv);
% 
% m = 1000;
% f1 = 0;
% f2 = 20;
% w = exp(-j*2*pi*(f2-f1)/(m*Fs));
% a = exp(j*2*pi*f1/Fs);
% z = czt(h,m,w,a) ./ (m * Fs/100);
% 
% fnz = (0:m-1)'/m;
% fz = (f2-f1)*fnz + f1;
% 
% figure (1) ; clf
% subplot(2,1,1) ; hold on
%     plot(fz,2*abs(z),'r')
%     plot(fy,2*abs(y),'b')
%     xlim([f1 f2])
%     legend('FFT','CZT')
% subplot(2,1,2) ; hold on
%     plot(fz,angle(z),'r')
%     plot(fy,angle(y),'b')
%     xlim([f1 f2])
%     xlabel('Frequency (Hz)')

end