function [daq_signal] = servo_daq(daq_fs, servo_signal, servo_time, servo_delay_range, servo_update_rate)
%% servo_signal: 
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

% servo_signal = replay.pos.body(:,2);
% servo_time = replay.time;

% servo_signal = 180*ones(size(servo_signal))

% Check inputs
if nargin < 5
    servo_update_rate = 50; % [hz]
    if nargin < 4
        servo_delay_range = [1 2]; % [ms]
    end
end

% Check servo signal
if range(servo_signal) > 180 % can't move more than 180 deg
    error('servo can only move within 180 deg')
end

if any(servo_signal < 0) % can't go below 0 deg
    servo_signal = servo_signal - min(servo_signal);
end

% Map servo signal to PWM pulse times
servo_pulse_signal = interp1([0 180], servo_delay_range, servo_signal);

% Map servo signal to servo update rate
servo_ts = 1/ servo_update_rate;
servo_update_time = (0:servo_ts:servo_time(end))';
servo_update_signal = interp1(servo_time, servo_pulse_signal, servo_update_time);

daq_ts = 1/ daq_fs;
daq_time = (0:daq_ts:servo_time(end))';
daq_signal = nan(size(daq_time));

% Generate PWM control signal
n_point = length(servo_update_signal);
for n = 1:n_point-1
    if n == 1
        start_time = servo_update_time(n);
        [~, daq_startI] = min(abs(start_time - daq_time));
    else
        daq_startI = daq_endI + 1;
    end
 	end_time = servo_update_time(n+1);
    [~, daq_endI] = min(abs(end_time - daq_time));
    
    pulse_length = round(daq_fs * servo_update_signal(n) * 1e-3);
    pulse_end = daq_startI + pulse_length;
    daq_signal(daq_startI:pulse_end) = 1;
    daq_signal((pulse_end + 1):daq_endI) = 0;
end

dcheck = diff(daq_signal);
dcheck = [dcheck ; dcheck(end)];
[~,pks] = findpeaks(abs(dcheck));

dpks = diff(pks);
max_d = daq_fs * 2e-3
dpks(dpks > max_d) = [];

end