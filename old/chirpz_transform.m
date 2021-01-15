function [Z,fz,mag,phs] = chirpz_transform(x,t,f1,f2,nfft)
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
if length(t) == 1
    fs = t;
else
    fs = 1/(mean(diff(t)));
end

% Number of points in DFT
if (nargin < 5) || isempty(nfft)
    T = t(end);
    nfft = ceil((f2 - f1)/(1/T));
end
m = nfft + 1;

% Make frequency vector
fnz = (0:m-1)'/m;
fz = (f2-f1)*fnz + f1;

% Compute Chirp-Z Transform
w = exp(-1j*2*pi*(f2-f1)/(nfft*fs));
a = exp(1j*2*pi*f1/fs);
Z = czt(x,m,w,a) ./ (m * fs/100);

% Get magniyude & phase
mag = abs(Z);
phs = angle(Z);

end