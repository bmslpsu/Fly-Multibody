function [Pxx,F,RPxx,Fc] = computeperiodogram2(x,win,nfft,esttype,Fs,options)
%COMPUTEPERIODOGRAM   Periodogram spectral estimation.
%   This function is used to calculate the Power Spectrum Sxx, and the
%   Cross Power Spectrum Sxy.
%
%   Sxx = COMPUTEPERIODOGRAM(X,WIN,NFFT) where x is a vector returns the
%   Power Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Sxy = COMPUTEPERIODOGRAM({X,Y},WIN,NFFT) returns the Cross Power
%   Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Inputs:
%    X           - Signal vector or a cell array of two elements containing
%                  two signal vectors.
%    WIN         - Window
%    NFFT        - Number of frequency points (FFT) or vector of
%                  frequencies at which periodogram is desired
%    ESTTYPE     - A string indicating the type of window compensation to
%                  be done. The choices are:
%                  'ms'    - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'power' - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'psd'   - compensate for Power Spectral Density (PSD);
%                            maintain correct area under the PSD curve.
%     REASSIGN   - A logical (boolean) indicating whether or not to perform
%                  frequency reassignment
%
%   Output:
%    Sxx         - Power spectrum [Power] over the whole Nyquist interval.
%      or
%    Sxy         - Cross power spectrum [Power] over the whole Nyquist
%                  interval.
%
%    F           - (vector) list frequencies analyzed
%    RSxx        - reassigned power spectrum [Power] over Nyquist interval
%                  has same size as Sxx.  Empty when 'reassigned' option
%                  not present.
%    Fc          - center of gravity frequency estimates.  Same size as
%                  Sxx.  Empty when 'reassigned' option not present.
%

%   Copyright 1988-2019 The MathWorks, Inc.
%#codegen

% Copied from signal/private/computeperiodogram with some modifications:
% 1. Allows ND-array computations
% 2. Doesn't include the branches for two input signals or reassigned
% periodograms.

narginchk(5,7);
if nargin>5
    reassign = options.reassign;
    assert(~reassign, "Reassigned periodogram computation is not available.");
end

% use normalized frequencies when Fs is empty
if isempty(Fs) || isnan(Fs)
    Fs = 2*pi;
end

% Validate inputs and convert row vectors to column vectors.
[x1,~,~,~,win1] = validateinputs(x,win,nfft);

% Window the data
xw = x1.*win1;

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels

% [Xx,F] = computeDFT2(xw,nfft,Fs);
% t = (0:(1/Fs):(length(x1)-1)*(1/Fs))';
% T = t(end);
f1 = 0;
f2 = Fs/2;
% nfft = ceil((f2 - f1)/(1/T));
[Xx,F] = chirpz(xw, Fs, f1, f2, length(nfft)-1);

% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if any(strcmpi(esttype,{'ms','power'}))
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win1)^2;
else
    U = win1'*win1;  % compensates for the power of the window.
end

% We do a real cast to make Pxx real in codegen.
% Pxx = real(Xx.*conj(Xx)/U);                % Auto spectrum.
Pxx = Xx;

RPxx = cast([],'like',Pxx);
Fc = [];
end



%--------------------------------------------------------------------------
function [x1,Lx,y,is2sig,win1] = validateinputs(x,win,~)
% Validate the inputs to computexperiodogram and convert row vectors to
% column vectors for backwards compatiblility with R2014a and prior
% releases

% Set defaults and convert to row vectors to columns.
is2sig= false;
win1   = win(:);
Lw    = length(win);

% Two signal vectors are not supported.
assert(~iscell(x), "Two signal inputs are not supported.");
y = [];

validateattributes(x,{'single','double'}, {'finite','nonnan'},'periodogram','x');
if isvector(x)
    x1 = x(:);
else
    x1 = x;
end

Lx = size(x1,1);

% ND-arrays are supported here.
coder.internal.errorIf(Lx ~= Lw, 'signal:computeperiodogram:invalidWindow', 'window');

end

% -------------------------------------------------------------------------

% LocalWords:  Sxx Sxy NFFT ESTTYPE RSxx Fc Fs Yy computexperiodogram
% LocalWords:  compatiblility nonnan
