function [Z,fz,mag,phs] = chirpz(x,t,f1,f2,nfft)
%% chirpz: Computes DFT using Chirp Z-Transform
%
%   INPUT:
%       x       : 	data
%       t       :  	time vector [s] or sampling frequency [Hz]
%       f1      :   lower frequency bound
%       f2      :   upper frequency bound
%       nfft   :   # of frequency bins
%
%   OUTPUT:
%       Z   :   Chirp Z-Transform
%       Fz  :  	frequency vector [Hz]
%       Mag : 	magnitude vector
%       Phs : 	phase vector [rad]
%

% Sampling frequency [Hz]
if length(t) == 1
    fs = t;
else
    fs = 1/(mean(diff(t)));
end
t = (0:(1/fs):(length(x)-1)*(1/fs))';
T = t(end);
fv = (0:(1/(1*T)):fs/2)';
n = length(fv);

% Number of points in DFT
if (nargin < 5) || isempty(nfft)
    nfft = ceil((f2 - f1)/(1/T));
else
    n = nfft;
end
m = nfft + 1;

% Make frequency vector
fnz = (0:m-1)'/nfft;
fz = (f2-f1)*fnz + f1;

% Compute Chirp-Z Transform
w = exp(-1j*2*pi*(f2-f1)/(nfft*fs));
a = exp(1j*2*pi*f1/fs);
Z = czt(x,m,w,a) ./ n;

% Get magniyude & phase
mag = abs(Z);
phs = angle(Z);

end