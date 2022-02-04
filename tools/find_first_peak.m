function [first_peak] = find_first_peak(signal, time, Fc, thresh, showplot)
%% find_first_peak: get the start of the first peak in signal that meets criteria
%
%   INPUT:
%     	signal        	: signal vector
%     	time          	: time vector
%       Fc              : cut off frequency [hz]
%       thresh       	: velocity threshold
%       showplot        : show plot (boolean)
%       
%   OUTPUT:
%       first_peak          : signal to output from DAQ
%

if nargin < 5
    showplot = true;
end

Fs = 1 / mean(diff(time));

% Filter & take derivative of signal
if ~isempty(Fc)
    [b, a] = butter(5, Fc / (2*Fs), 'low');
    signal_filt = filtfilt(b, a, signal);
    signal_vel = central_diff(signal_filt, 1/Fs);
    signal_vel = filtfilt(b, a, signal_vel);
else
    signal_filt = signal;
    signal_vel = central_diff(signal_filt, 1/Fs);
end

% Find first peak location
% [~, locs] = findpeaks(signal_vel, 'MinPeakHeight', thresh);
% first_peak = locs(1);
first_peak = find(signal_vel > thresh, 1, 'first');

% Find where efirst peak starts
start_thresh = 0.10*signal_vel(first_peak);
search_signal = signal_vel(1:first_peak);
search_time = time(1:first_peak);
start_index = find(search_signal < start_thresh, 1, 'last');

if showplot
    subplot(2,1,1) ; cla ; hold on
        plot(time, signal, 'k', 'LineWidth', 1)
        plot(time(first_peak), signal(first_peak), '.r', 'MarkerSize', 20)
        plot(search_time, signal(1:first_peak), 'b', 'LineWidth', 1)
        plot(time(start_index), signal(start_index), '.g', 'MarkerSize', 20)
    subplot(2,1,2) ; cla ; hold on
        plot(time, signal_vel, 'k', 'LineWidth', 1)
        plot(search_time, search_signal, 'b', 'LineWidth', 1)
        plot(time(first_peak), signal_vel(first_peak), '.r', 'MarkerSize', 20)
        plot(time(start_index), signal_vel(start_index), '.g', 'MarkerSize', 20)
        yline(thresh, '--m');
end

end