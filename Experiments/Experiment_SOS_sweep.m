function [] = Experiment_SOS_sweep(Fn)
%% Experiment_SOS_norm: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_SOS_vel_v1';

%% EXPERIMENTAL PARAMETERS
n_tracktime = 21;           % length(func)/fps; seconds for each EXPERIMENT
% n_resttime = 1;         	% seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 1;                  % pattern ID
yPos = 5;                   % 30 deg spatial frequency
% funcID = 1;             	% position function ID
xUpdate = 50;               % function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 0:2;                	% Analog input channels
AO = 1;                     % Analog output channels

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
[s,~] = MC_USB_1208FS_PLUS(Fs,AI,AO);

%% Camera Trigger Signal
off = 0.1;
t = 0:1/s.Rate:n_tracktime + off;
TRIG = ((1/2)*(square(2*pi*FPS*t,5) - 1)');
TRIG(TRIG==-1) = 4;
% start_off = round(Fs * 0.5*(1/FPS));
end_off = round(Fs*off);
% TRIG(1:start_off) = 0;
TRIG(end-end_off:end) = 0;

% subplot(1,2,1)
% plot(t,TRIG)
% ylim([-0.1 4.1])
% xlim([-0.01 0.1])
% 
% subplot(1,2,2)
% plot(t,TRIG)
% ylim([-0.1 4.1])
% xlim([t(end) - 0.1 , t(end)])

%% Camera Setup
[vid,src] = Basler_acA640_750um(nFrame);

%% Set variable to control positionn function
func = [1 2 3]';  % position function indicies
vel = [62 103 148]'; % normalizee veloctiy of each SOS function
n_func = length(func);	% # of functions

% Create sequence of randomly shuffled frequencies
func_all = nan(n_func*n_rep,1);
pp = 0;
for kk = 1:n_rep
    func_rand = func(randperm(n_func),:);    % reshuffle randomly
    func_all(pp+1:pp+n_func,1) = func_rand;  % add rep
    pp = kk*n_func;
end

n_trial = n_rep * n_func;

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   Vel = %i \n', ii, vel(func_all(ii)))
    preview(vid) % open video preview window
    
    Panel_com('stop')
    
    % Set AO trigger to 0
 	queueOutputData(s,zeros(5000,1)) % set trigger AO signal to 0
    [~,~] = s.startForeground; % data collection
    
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
	
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    queueOutputData(s,TRIG) % set trigger AO signal
    T = timer('StartDelay',0.5,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        stop(vid) % stop video buffer
        Panel_com('stop') % stop stimulus
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % SPIN BUFFER
    Arena_Ramp(1,10)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_vel_' ...
        num2str(vel(func_all(ii))) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end