function [pulse_dir_signal, pulse_signal, dir_signal, daq_time] = stepper_signal_daq(daq_fs, signal, time, step_size, showplot)
%% stepper_signal_daq: generates analog output signal for DAQ, to send to arduino to control stepper motor
%
%   INPUT:
%     	daq_fs              : DAQ sampling frequency [hz]
%       signal              : desired stepper trajectory [deg]
%       time                : desired stepper time [s]
%       step_size           : step size [deg] (adjust for micro stepping)
%       showplot            : show plot (boolean)
%       
%   OUTPUT:
%       pulse_dir_signal	: signal that encodes pulses & direction
%       pulse_signal     	: encodes pulses
%       daq_time            : encodes direction
%       daq_time            : DAQ time
%

if nargin < 5
    showplot = true;
end

medV = 2.5;

% Convert to DAQ sampling rate time
daq_ts = 1/ daq_fs;
daq_time = (0:daq_ts:time(end))';
daq_signal = interp1(time, signal, daq_time);
daq_signal = daq_signal - daq_signal(1);

% Round signal to nearest step size
step_signal = step_size*round(daq_signal./step_size);

diff_signal = diff(step_signal);
diff_signal = [diff_signal ; diff_signal(end)];

% Make sure every step can be set by the DAQ
check_range = unique(round(diff_signal,5));
assert(length(check_range) <= 3, 'DAQ sampling rate to low to achieve desired trajectory')

% Find each rising edge in either direction
[~, locs] = findpeaks(abs(diff_signal), 'MinPeakHeight', mean(abs(diff_signal)));

% Extend the peak width
n_peak = length(locs);
diff_signal_wide = diff_signal;
tol = 2;
for n = 1:n_peak-1
    k1 = locs(n);
    k2 = locs(n + 1);
    dk = k2 - k1;
    if dk > tol % more than 2 data points between rising edges
        win = dk - tol; % extend spike by this many points
        val = diff_signal(locs(n)); % use this value
        diff_signal_wide(k1:(k1+win)) = val; % modify signal
    else
        % do nothing
    end
end

% Map signals to 0-5V
pulse_dir_signal = medV*diff_signal_wide./max(diff_signal_wide) + medV; % encode puleses & direction
pulse_signal = abs(diff_signal_wide); % just encode pulses
pulse_signal = 2*medV*pulse_signal./max(pulse_signal);

% Make direction signal
dir_signal = nan(size(pulse_signal));
dir_signal(1:locs(1)) = sign(diff_signal(locs(1)));
for n = 1:n_peak-1
    k1 = locs(n);
    k2 = locs(n + 1);
    dk = k2 - k1;
    mid = floor(dk/2);
    pk1 = sign(diff_signal(k1));
    pk2 = sign(diff_signal(k2));
    if pk1 == pk2 % pulses are same direction
        dir_signal(k1:k2) = pk1;
    else % pulses are opposite direction
        dir_signal(k1:(k1+mid)) = pk1;
        dir_signal((k1+mid)+1:k2) = pk2;
    end
end
dir_signal(locs(end):end) = sign(diff_signal(locs(end)));
dir_signal = medV * (dir_signal./max(dir_signal)) + medV;

if showplot
    close all
    ax(1) = subplot(3,1,1) ; cla ; hold on
        plot(daq_time, daq_signal, 'k', 'LineWidth', 1)
        plot(daq_time, step_signal, 'r', 'LineWidth', 1)
    
    ax(2) = subplot(3,1,2) ; cla ; hold on
        plot(daq_time, diff_signal, 'k', 'LineWidth', 0.5)
        plot(daq_time, diff_signal_wide, 'r', 'LineWidth', 0.5)
    
    ax(3) = subplot(3,1,3) ; cla ; hold on
        plot(daq_time, pulse_dir_signal, 'b', 'LineWidth', 0.5)
        %plot(daq_time, pulse_signal, 'g', 'LineWidth', 0.5)
        plot(daq_time, 1.2*dir_signal - 0.1, 'r', 'LineWidth', 0.5)
        
    ylim([-0.2 5.2])
    yticks([0 2.5 5])
    linkaxes(ax, 'x')
end

end