function [pkFreq,pkMag,pkPhase,pkIdx] = getfreqpeaks(freq,mag,phase,IOFreq,fTol,debug)
%% getfreqpeaks: get magnitude & phase values of frequency domain data at input-output (IO) freqeuncies
%
%   INPUT:
%       freq        : frequency vector
%       mag         : magnitude vector
%       phase       : Phase vector
%       IOFreq     	: unique freuncies present (if empty, find peaks automatically)
%       fTol        : frequency search tolerance (default=2*T) >>> varargin{1}
%       debug       : showplot (default=false) >>> varargin{2}
%
%   OUTPUT:
%       pkFreq      : IO frequencies
%       pkMag       : manitude at uFreq
%       pkPhase 	: phase at uFreq
%       pkIdx       : index of peaks
%

if isempty(fTol)
    fTol = 2*mean(diff(freq)); % default tolerance range
end

if isempty(phase)
    phase = zeros(size(freq));
end

if nargin<6
    debug = false; % default is off
end

if isempty(IOFreq)
    autoflag = true; % IO frequencies not provided
    
    % Remove noise based on threshold
    medMag = median(mag);
    stdMag = std(mag);
    thresh = medMag + 2*stdMag;
    magtest = mag;
    magtest(mag<thresh) = 0;
    
    % Detetect peaks
	[~,locs,~,~] = findpeaks(magtest);
    IOFreq = freq(locs);
else
    autoflag = false; % IO freuencies provided
end

nn          = length(freq);     % of # of frequency domain points
nFreq       = length(IOFreq);  	% # of frequencies to use
Ff          = nan(nFreq,1);     % closest frequency that max mag occurs
pkMag       = nan(nFreq,1);     % magnitude at frequencies
pkPhase     = nan(nFreq,1);     % phase at frequencies
pkIdx     	= nan(nn,nFreq);    % indicies at frequencies

for kk = 1:nFreq
    fRange = [IOFreq(kk)-fTol , IOFreq(kk)+fTol]; % frequency search range
    idxRange = freq>=fRange(1) & freq<=fRange(2); % index search range
    magRange = mag(idxRange); % magnitude search range
    fIdx = mag==max(magRange); % index of max magnitude
    Ff(kk) = freq(fIdx);
    
    pkMag(kk) = max(magRange); % store magnitude at uFreq
    pkPhase(kk) = phase(fIdx); % store phase at uFreq
    
    pkIdx(:,kk) = fIdx;
end

pkIdx = logical(sum(pkIdx,2));
if sum(pkIdx,1)<nFreq
    error('Error: conflicting frequencies')
end
[pkIdx,~] = find(pkIdx==true);
pkFreq = freq(pkIdx);

% Debug plot
if debug % plot magnitude & phase
fig = figure ; clf
fig.Color = 'w';
    subplot(2,1,1) ; hold on ; ylabel('Magnitude')
    plot(freq,mag,'*-b')
    if autoflag; plot(freq,magtest,'y') ; end
    plot(Ff,pkMag,'or')
    xlim([0 max(IOFreq+1)])

    subplot(2,1,2) ; hold on ; ylabel('Phase')
    plot(freq,phase,'*-b')
    plot(Ff,pkPhase,'or')
    xlim([0 max(IOFreq+1)])
    xlabel('Frequency')
end
end