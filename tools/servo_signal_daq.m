function [daq_signal, daq_time] = servo_signal_daq(daq_fs, servo_signal, servo_time, showplot)
%% servo_signal_daq: 
%
%   INPUT:
%     	daq_fs           	: DAQ sampling frequency [hz]
%       servo_signal        : desired servo trajectory [deg]
%       servo_time          : desired servo time [s]
%       servo_delay_range  	: [min max] delay range for mapping angles to PWM pulse width [ms]
%       servo_update_rate  	: how often to send PWM pulses (default is 50 hz)
%       
%   OUTPUT:
%       daq_signal       	: PWM signal to output from DAQ
%

if nargin < 4
    showplot = true;
end

% Check servo signal
if range(servo_signal) > 180 % can't move more than 180 deg
    error('servo can only move within 180 deg')
end

if any(servo_signal < 0) % can't go below 0 deg
    servo_signal = servo_signal - min(servo_signal);
end

% Map servo signal to PWM pulse times
deg_range = [min(servo_signal) max(servo_signal)];
servo_pulse_signal = interp1(deg_range, [0 5], servo_signal);

daq_ts = 1/ daq_fs;
daq_time = (0:daq_ts:servo_time(end))';
daq_signal = interp1(servo_time, servo_pulse_signal, daq_time);

if showplot
   plot(daq_time, daq_signal, 'k', 'LineWidth', 1)
end

end