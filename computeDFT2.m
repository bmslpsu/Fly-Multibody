function [Xx,f] = computeDFT2(xin,nfft,varargin)
%#codegen
%COMPUTEDFT Computes DFT using FFT or Goertzel
%   This function is used to calculate the DFT of a signal using the FFT
%   or the Goertzel algorithm.
%
%   [XX,F] = COMPUTEDFT(XIN,NFFT) where NFFT is a scalar and computes the
%   DFT XX using FFT. F is the frequency points at which the XX is
%   computed and is of length NFFT.
%
%   [XX,F] = COMPUTEDFT(XIN,F) where F is a vector with at least two
%   elements computes the DFT XX using the Goertzel algorithm.
%
%   [XX,F] = COMPUTEDFT(...,Fs) returns the frequency vector F (in hz)
%   where Fs is the sampling frequency
%
%   Inputs:
%   XIN is the input signal
%   NFFT if a scalar corresponds to the number of FFT points used to
%   calculate the DFT using FFT.
%   NFFT if a vector corresponds to the frequency points at which the DFT
%   is calculated using goertzel.
%   FS is the sampling frequency

% Copyright 2006-2019 The MathWorks, Inc.

% [1] Oppenheim, A.V., and R.W. Schafer, Discrete-Time Signal Processing,
% Prentice-Hall, Englewood Cliffs, NJ, 1989, pp. 713-718.
% [2] Mitra, S. K., Digital Signal Processing. A Computer-Based Approach.
% 2nd Ed. McGraw-Hill, N.Y., 2001.

% Copied from signal/private/computeDFT with the required modifications to:
% 1. Allow ND-array computations with the Goertzel algorithm.


narginchk(2,3);
if nargin > 2
    Fs = varargin{1};
else
    Fs = 2*pi;
end

nx = size(xin,1);

isfreqScalar = isscalar(nfft);

if isfreqScalar
    [Xx,f] = computeDFTviaFFT(xin,nx,nfft(1),Fs);
else
    f = nfft(:); % if nfft is a vector then it contains a list of freqs
    
    % see if we can get a uniform spacing of the freq vector
    [fstart, fstop, m, maxerr] = signal.internal.spectral.getUniformApprox(f);
    
    % see if the ratio of the maximum absolute deviation relative to the
    % largest absolute in the frequency vector is less than a few eps
    isuniform = maxerr < 3*eps(class(f));
    
    % check if the number of steps in Goertzel ~  1 k1 N*M is greater
    % than the expected number of steps in CZT ~ 20 k2 N*log2(N+M-1)
    % where k2/k1 is empirically found to be ~80.
    n = size(xin,1);
    islarge = m > 80*log2(nextpow2(m+n-1));
    
    if isuniform && islarge
        % use CZT if uniformly spaced with a significant number of bins
        Xx = computeDFTviaCZT(xin,fstart,fstop,m,Fs);
    else
        % use Goertzel if small number of bins or not uniformly spaced
        Xx = computeDFTviaGoertzel(xin,f,Fs);
    end
end

end

%-------------------------------------------------------------------------
function [Xx,f] = computeDFTviaFFT(xin,nx,nfft,Fs)
% Use FFT to compute raw STFT and return the F vector.

% Handle the case where NFFT is less than the segment length, i.e., "wrap"
% the data as appropriate.
xin_ncol = size(xin,2);
xin_nchan = size(xin,3);
xw = zeros(nfft,xin_ncol,xin_nchan,'like',xin);
if nx > nfft
    for j = 1:xin_ncol*xin_nchan
        wrappedData = datawrap(xin(:,j),nfft);
        xw(:,j) = wrappedData(:); %coder size inference:- infer it as Column Vector to match LHS
    end
else
    xw = xin;
end

Xx = fft(xw,nfft);
f = psdfreqvec('npts',nfft,'Fs',Fs);

end

%--------------------------------------------------------------------------
function Xx = computeDFTviaGoertzel(xin,f,Fs)
% Use Goertzel to compute raw DFT

f = mod(f,Fs);    % 0 <= f < = Fs
xm = size(xin,1); % NFFT

isInMATLAB = coder.target('MATLAB');
isInputComplex = ~isreal(xin);
isInputSingle = isa(xin,'single');

% wavenumber in cycles/period used by the Goertzel function
% (see equation 11.1 pg. 755 of [2])
if isInputSingle
    k = single(f/Fs*xm);
else
    k = f/Fs*xm;
end

% Goertzel algorithm is defined per segment/window:
sizeIn = size(xin);
Xx = zeros([numel(f) sizeIn(2:end)], 'like', xin);
for ii = 1:sizeIn(2)
    out = signal.internal.goertzel.callGoertzel(squeeze(xin(:, ii, :)),k,isInputComplex,isInputSingle,isInMATLAB);
    % Reshape back to nfft-by-1-by-numChannels
    Xx(:, ii, :) = permute(out, [1 3 2]);
end

end

%--------------------------------------------------------------------------
function Xx = computeDFTviaCZT(xin,fstart,fstop,npts,Fs)
% Use CZT to compute raw DFT

% start with initial complex weight
Winit = exp(2i*pi*fstart/Fs);

% compute the relative complex weight
Wdelta = exp(2i*pi*(fstart-fstop)/((npts-1)*Fs));

% feed complex weights into chirp-z transform
Xx = czt(xin, npts, Wdelta, Winit);

end

% LocalWords:  DFT XIN NFFT Fs hz FS Oppenheim Schafer Englewood Mitra nd Graw
% LocalWords:  nfft STFT npts wavenumber
