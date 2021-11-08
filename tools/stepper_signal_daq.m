function [daq_signal, daq_time] = stepper_signal_daq(daq_fs, signal, time, step_size, showplot)
%% stepper_signal_daq: generates analog output signal for DAQ, to send to arduino to control stepper motor
%
%   INPUT:
%     	daq_fs        	: DAQ sampling frequency [hz]
%       signal          : desired stepper trajectory [deg]
%       time            : desired stepper time [s]
%       step_size    	: step size [deg] (adjust for micro stepping)
%       showplot        : show plot (boolean)
%       
%   OUTPUT:
%       daq_signal    	: signal to output from DAQ
%       daq_time    	: DAQ time
%

if nargin < 5
    showplot = true;
end

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
assert(length(check_range) == 3, 'DAQ sampling rate to low to achieve desired trajectory')

% Map signal to 0-5V
diff_signal = 2.5*diff_signal./max(diff_signal) + 2.5;

if showplot
    close all
    subplot(2,1,1) ; cla ; hold on
    plot(daq_time, daq_signal, 'k', 'LineWidth', 1)
    plot(daq_time, step_signal, 'r', 'LineWidth', 1)
    
    subplot(2,1,2) ; cla ; hold on
    plot(daq_time, diff_signal, 'b', 'LineWidth', 0.5)
    ylim([-0.2 5.2])
    yticks([0 2.5 5])
end

end