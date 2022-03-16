function [] = Experiment_SS_sweep_stepper_motor_drum(Fn)
%% Experiment_SS_sweep_stepper_motor_drum: runs a experiment using the LED arena and fly panel
% Fn is the fly number
clc
daqreset
imaqreset
% Fn = 0;

%% Set directories & experimental parameters
% root = 'C:\BC\Experiment_SS_vel_250_motor_drum';
root = 'C:\BC\Experiment_SS_vel_250_motor_passive';

step_size = 1.8/16;

%% EXPERIMENTAL PARAMETERS
n_tracktime = 10 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
% n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
% patID = 1;                  % pattern ID
% xUpdate = 400;           	% function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 10000;              	% DAQ sampling rate [Hz]

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
[s,~] = MC_USB_1208FS_PLUS_motor_all(Fs);
% [s.camera,~] = NI_USB_6212_motor(Fs,[],0);

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

%% Custom
% Replay
replay_path = fullfile(root, 'replay', 'replay_SS_vel_250_position.mat');
replay = load(replay_path);
replay_time = replay.replay.time;
replay_ts = mean(diff(replay_time));
replay_fs = 1 / replay_ts;
replay_time_new = (0:replay_ts:total_time)';

replay_body = replay.replay.pos.body_sine;
replay_body = [replay_body ; replay_body(1:replay_fs*0,:)];

replay_body(:,end) = 2*replay_body(:,end);
n_replay = size(replay_body,2);

start_pad = replay_body(1)*ones(round(0.5*replay_fs),n_replay);
replay_body_new = [start_pad ; replay_body];

current_size = size(replay_body_new,1);
new_size = length(replay_time_new);
size_increase = new_size - current_size + 1;
replay_body_new(current_size:new_size,:) = repmat(replay_body_new(end,:), [size_increase,1]);
motor_signal = replay_body_new;
motor_time = replay_time_new;

% close
% plot(replay_time_new, replay_body_new)

% clear pulse_signal dir_signal
for n = 1:n_replay
    [~, pulse_signal(:,n), dir_signal(:,n), ~] = stepper_signal_daq(s.Rate, ...
        motor_signal(:,n), motor_time, step_size, false);
end

%% Set variable to control position function
val = [replay.replay.IOFreq{:}]'; % amplitude of each SOS function in order in PControl
name = 'freq'; % name of identifier at end of file name

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

%% EXPERIMENT LOOP
clc
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   Freq = %f \n', ii, val_all(ii))
    preview(vid) % open video preview window
        
    % Set AO trigger to 0
    queueOutputData(s,zeros(s.Rate,3))
    [data, t_p] = s.startForeground;
    
    pause(1) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    AO = [pulse_signal(:, func_all(ii)), dir_signal(:, func_all(ii)), TRIG]; % set stepper motor control signal
    
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    queueOutputData(s, AO) % set stepper AO signal
	
    %T = timer('StartDelay', 0.5, 'TimerFcn', @(src,evt) Panel_com('start'));
    %start(T)
    
 	[data, t_p] = s.startForeground; % data collection
    stop(vid) % stop video buffer
    [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_' name '_' ...
        num2str(val_all(ii)) '.mat'];
    save(fullfile(root,fname), '-v7.3', 'data', 't_p', 'vidData', 't_v',...
        'step_size', 'AO', 'motor_signal', 'motor_time');
end

delete(vid)
disp('Done');
daqreset
imaqreset
end