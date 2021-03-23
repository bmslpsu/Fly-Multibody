close all
clear
clc

fs = 100;
fn = fs/2;	% nyquist frequency [Hz]
T = 20;
t = (0:(1/fs):T)';
n = length(t);
a = 10;
f = 1;
x = 0.1*t.^(2) + a*cos(2*pi*f*t) + 0*cos(4*pi*f*t);

window = round(10*fs);
percent_overlap = 0.8;
noverlap = round(percent_overlap*window);
fv = (0:(1/(1*T)):fn)';
nfft = length(fv)-1;
% [welchOut] = welch(x, 'pwelch', window, noverlap, [], fs);
% pxx = pwelch(x);
funcName = 'pwelch';
[welchOut,fw] = welch2(x, funcName, window, noverlap, fv, fs, 'power');

% [Xx,f] = computeDFT2(x, nfft, fs);
% [~, ~ , ~ , FREQ] = FFT(t,x);

f1 = 0;
f2 = fs/2;
[Z,fz] = chirpz(x,fs,f1,f2);

% hold on
% plot(fz, abs(Z))
% plot(fw, 7.42*abs(welchOut))
% xlim([0 5])

fig = figure (1);
set(fig, 'Color', 'w')
clear ax
ax(1) = subplot(3,1,1); hold on
    plot(t, x, 'b')
ax(2) = subplot(3,1,2); hold on
    plot(fz, abs(Z), 'k')
    plot(fw, 7.42*abs(welchOut)./2, 'r')
ax(3) = subplot(3,1,3); hold on
    plot(fz, angle(Z), 'k')
    plot(fw, angle(welchOut), 'r')

set(ax, 'LineWidth', 1)
set(ax(2:3), 'XLim', [0 5])

%%
% fv = (linspace(0, 1, fix(n/2)+1)*fn)';  % frequency vector [Hz]
fts = welchOut;                       % normalized fourier transform
Iv = 1:length(fv);                  	% index vector
mag = abs(fts(Iv));                   % magnitude
phs = angle(fts(Iv));                   % phase [rad]

subplot(2,1,1) ; hold on
    plot(fv,mag)
subplot(2,1,2) ; hold on
    plot(fv,phs)

%%
close all
clear
clc

fs = 50;
fn = fs/2;	% nyquist frequency [Hz]
T = 10;
t = (0:(1/fs):T)';
n = length(t);
a = 10;
f = 1;
x = 2*t + a*cos(2*pi*f*t) + 0*cos(4*pi*f*t);
win = hann(n);
% win = hamming(n);
y = win.*x;

f1 = 0;
f2 = fs/2;
[X,fx] = chirpz(x, fs, f1, f2);
[Y,fy] = chirpz(y, fs, f1, f2);

fig = figure (1);
set(fig, 'Color', 'w')
clear ax
ax(1) = subplot(3,1,1); hold on
    plot(t, x, 'k')
    plot(t, y, 'r')
ax(2) = subplot(3,1,2); hold on
    plot(fx, abs(X), 'k')
    plot(fy, 2*abs(Y), 'r')
ax(3) = subplot(3,1,3); hold on
    plot(fx, angle(X), 'k')
    plot(fy, angle(Y), 'r')

set(ax, 'LineWidth', 1)

%%
% % fs = 1000;
% % d = designfilt('lowpassfir','FilterOrder',30,'CutoffFrequency',125, ...
% %     'DesignMethod','window','Window',@rectwin,'SampleRate',fs);
% % h = tf(d);
% fs = All.Fs;
% h = All.X_step;
% T = All.time(end);
% fn = fs/2;
% L = length(h);
% res = 1/10;
% 
% n = fn/res;
% y = fft(h,n);
% 
% fy = (linspace(0, 1, fix(L/2)+1)*fn)';  % frequency vector [Hz]
% fts = fft(h)/L;                        	% normalized fourier transform
% Iv = 1:length(fy);                  	% index vector
% mag = 2*abs(fts(Iv));                   % magnitude
% phs = angle(fts(Iv));                   % phase [rad]
% FREQ = fts(Iv);                         % complex frequency domain data
% 
% % m = 1000;
% f1 = 0;
% f2 = 50;
% m = 1*(f2 - f1)/res;
% w = exp(-1j*2*pi*(f2-f1)/(m*fs));
% a = exp(1j*2*pi*f1/fs);
% z = czt(h,m,w,a);
% 
% fn = (0:m-1)'/m;
% % fy = fs*fn;
% fz = (f2-f1)*fn + f1;
% 
% hold on
% plot(fy, mag, 'b')
% plot(fz, 2*abs(z)./L, 'r')
% xlim([0 10])
% legend('FFT','CZT')
xlabel('Frequency (Hz)')