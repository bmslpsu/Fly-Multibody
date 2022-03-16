function [] = Experiment_SOS_sweep_stepper_motor_visual(Fn)
%% Experiment_SOS_sweep_stepper_motor_visual: runs a experiment using the LED arena and fly panel
% Fn is the fly number
clc
daqreset
imaqreset
% Fn = 0;

%% Set directories & experimental parameters

yPos = 5;

switch yPos
    case 12
        root = 'C:\BC\Experiment_SOS_vel_v2_motor_bright';
    case 1
        root = 'C:\BC\Experiment_SOS_vel_v2_motor_dark';
    case 5
        root = 'C:\BC\Experiment_SOS_vel_v2_motor_grating';
    otherwise
        error('must pick other y-position')
end

val = [70]'; % amplitude of each SOS function in order in PControl
name = 'vel'; % name of identifier at end of file name

step_size = 1.8/16;

%% EXPERIMENTAL PARAMETERS
n_tracktime = 20 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 1;                  % pattern ID
xUpdate = 400;           	% function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 10000;              	% DAQ sampling rate [Hz]

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
[s,~] = MC_USB_1208FS_PLUS_motor_panel(Fs);

%% Camera Trigger Signal
off = 0.1;
total_time = n_tracktime + off;
t = 0:1/s.Rate:total_time;
TRIG = ((1/2)*(square(2*pi*FPS*t,5) - 1)');
TRIG(TRIG==-1) = 4;
end_off = round(Fs*off);
TRIG(end-end_off:end) = 0;

%% Camera Setup
[vid,src] = Basler_acA640_750um(nFrame);

%% Set variable to control position function
n_func = length(val); % # of functions
func = (1:n_func)'; % position function indicies

% Create sequence of randomly shuffled functions
func_all = nan(n_func*n_rep,1);
pp = 0;
for kk = 1:n_rep
    func_rand = func(randperm(n_func),:);    % reshuffle randomly
    func_all(pp+1:pp+n_func,1) = func_rand;  % add rep
    pp = kk*n_func;
end
val_all = val(func(func_all)); % true values (amplitude, velocity, etc.)
n_trial = n_rep * n_func;

% Replay
replay_path = fullfile(root, 'replay', 'replay_SOS_HeadFree_vel_v2_position_02-12-2021.mat');
replay = load(replay_path);
motor_pos = replay.replay.pos.body(:,2);
motor_time = replay.replay.time;

%% SOS
sos_path = fullfile(root, 'function', ...
    'ALL_position_function_SOS_Fs_400.63_T_20_vel_70_amp_1_1.54_2.38_3.67_5.66_8.72_13.45_20.75_32_freq_11.15_7.2_4.7_3.05_1.95_1.3_0.85_0.55_0.35.mat');
sos = load(sos_path);
% motor_pos = sos.All.X;
% motor_time = sos.All.time;

%% Visual input
ttt = (0:(1/s.Rate):total_time)';
test = 2 + 1.5*sin(2*pi*0.5*ttt);

%% Set stepper motor control signal
clc
pad_time = 0.5;

motor_ts = mean(diff(motor_time));
motor_fs = 1 / motor_ts;
motor_time_new = (0:motor_ts:total_time)';

d_off = motor_time_new(end) - total_time;
if d_off < 0
    motor_time_new(end+1) = total_time;
end

start_pad = motor_pos(1)*ones(round(pad_time*motor_fs),1);
motor_new = [start_pad ; motor_pos];
motor_new(end:length(motor_time_new)) = motor_new(end);

[~, pulse_signal, dir_signal, ~] = stepper_signal_daq(s.Rate, motor_new, motor_time_new, step_size, true);
AO = [pulse_signal , dir_signal, TRIG, test];

%% EXPERIMENT LOOP
clc
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   Amp = %i \n', ii, val_all(ii))
    preview(vid) % open video preview window
    
    Panel_com('stop')
    
    % Set AO trigger to 0
    queueOutputData(s,zeros(s.Rate,size(AO,2)))
    [data, t_p] = s.startForeground;
    
    pause(1) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    Panel_com('set_pattern_id', patID);	% set pattern
    pause(n_pause)
    Panel_com('set_position', [1, yPos]); % set starting position (xpos,ypos)
    pause(n_pause)
    Panel_com('set_posfunc_id',[1,func_all(ii)]); % arg1 = channel (x=1,y=2); arg2 = funcID
    pause(n_pause)
	Panel_com('set_funcX_freq', xUpdate); % update rate for x-channel
    pause(n_pause)
    Panel_com('set_funcY_freq', 50); % update rate for y-channel
    pause(n_pause)
    Panel_com('set_mode', [3,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    pause(n_pause)
    Panel_com('send_gain_bias', [0 0 0 0])
    Panel_com('start')
	
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    queueOutputData(s, AO) % set stepper AO signal
	
    %T = timer('StartDelay', 0.5, 'TimerFcn', @(src,evt) Panel_com('start'));
    %start(T)
    
 	[data, t_p] = s.startForeground; % data collection
    stop(vid) % stop video buffer
    Panel_com('stop') % stop stimulus
    [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    Panel_com('set_position', [1, 5])
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_' name '_' ...
        num2str(val_all(ii)) '.mat'];
    save(fullfile(root,fname), '-v7.3', 'data', 't_p', 'vidData', 't_v',...
        'step_size', 'AO', 'yPos', 'motor_new', 'motor_time_new');
end

delete(vid)
disp('Done');
daqreset
imaqreset
end