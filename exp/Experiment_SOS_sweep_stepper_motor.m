function [] = Experiment_SOS_sweep_stepper_motor(Fn)
%% Experiment_SOS_sweep_motor: runs a experiment using the LED arena and fly panel
% Fn is the fly number
clc
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_SOS_vel_v2_motor';
% val = [42 70 95]'; % amplitude of each SOS function in order in PControl
val = [70]'; % amplitude of each SOS function in order in PControl
name = 'vel'; % name of identifier at end of file name

step_size = 1.8/16;

%% EXPERIMENTAL PARAMETERS
n_tracktime = 20 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 1;                  % pattern ID
yPos = 5;                   % 30 deg spatial frequency
xUpdate = 400;           	% function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 5000;              	% DAQ sampling rate [Hz]
AI = 0:2;                	% Analog input channels
AO = 0:1;                 	% Analog output channels

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
[s,~] = MC_USB_1208FS_PLUS(Fs,AI,AO);

%% Camera Trigger Signal
off = 0.1;
t = 0:1/s.Rate:n_tracktime + off;
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

%% Set servo motor control signal
servo_time = (0:(1/1000):(n_tracktime + off))';
% servo_signal = 20*sin(2*pi*pi*1*servo_time);
servo_signal = 30*servo_time;
[motor_signal, daq_time] = stepper_signal_daq(Fs, servo_signal, servo_time, step_size, true);

% replay_path = fullfile(root, 'replay', 'replay_SOS_HeadFree_vel_v2_position_02-12-2021.mat');
% replay = load(replay_path);
% replay_body = replay.replay.pos.body(:,2);
% replay_time = replay.replay.time;
% 
% replay_ts = mean(diff(replay_time));
% replay_fs = 1 / replay_ts;
% replay_time_new = (0:replay_ts:(n_tracktime + off))';
% start_pad = replay_body(1)*ones(round(0.5*replay_fs),1);
% replay_body_new = [start_pad ; replay_body];
% replay_body_new(end:length(replay_time_new)) = replay_body_new(1);
% 
% [motor_signal, daq_time, max_range] = servo_signal_daq(Fs, replay_body_new, replay_time_new, true);

%% EXPERIMENT LOOP
clc
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   Amp = %i \n', ii, val_all(ii))
    preview(vid) % open video preview window
    
    Panel_com('stop')
    
    % Set AO trigger to 0
 	queueOutputData(s,zeros(5000,2))
    [~,~] = s.startForeground;
    
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
    Panel_com('set_mode', [4,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    pause(n_pause)
    Panel_com('send_gain_bias', [0 0 0 0])
	
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    AO = [motor_signal, TRIG];
    queueOutputData(s, AO) % set trigger AO signal
    %T = timer('StartDelay',0.5,'TimerFcn',@(src,evt) Panel_com('start'));
    %start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        stop(vid) % stop video buffer
        Panel_com('stop') % stop stimulus
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % SPIN BUFFER
    %Arena_Ramp(1,10)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_' name '_' ...
        num2str(val_all(ii)) '.mat'];
    %save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end