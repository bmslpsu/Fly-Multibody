function [signal] = steps2signal(pulse, dir, step_size, showplot)
%% steps2signal: converts sequence of steps & direction to time domain signal
%
%   INPUT:
%     	pulse        	: signal containing pulse information
%     	dir             : signal containing direction information
%       step_size    	: step size [deg] (adjust for micro stepping)
%       showplot        : show plot (boolean)
%       
%   OUTPUT:
%       signal          : signal to output from DAQ
%

if nargin < 4
    showplot = true;
    if nargin < 3
        step_size = [];
    end
end

% Convert pulse to logical & dir to sign
pulse = round(pulse);
pulse = pulse / max(pulse);
dir = round(dir);
dir = dir / (max(dir)/2);
dir = dir - 1;

% Make index & signal vectors
I = (1:length(pulse))';
signal = nan(size(I));

% Find pulse locations
[~, locs] = findpeaks(pulse, 'MinPeakHeight', 0.5);

% Reconstruct signal from pulses & direction
n_peak = length(locs);
signal(1:locs(1)-1) = 0;
for n = 2:n_peak
    signal(locs(n-1):locs(n)-1) = signal(locs(n-1)-1) + dir(locs(n))*1;
end
signal(locs(n_peak):end) = signal(locs(n_peak)-1);

% Convert to deg, otherwise leave in steps
if ~isempty(step_size)
    signal = signal * step_size;
end

if showplot
    plot(I, signal, 'k', 'LineWidth', 1)
    axis tight
end

end