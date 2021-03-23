function [] = Experiment_SS_sweep_rigid(Fn)
%% Experiment_SS_sweep_rigid: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_SS_vel_250_body_fixed_CL';
val = [3.75 7.5 11.25 18.75 26.25 41.25 60]'; % amplitude of each SS function in order in PControl
freq = [10.6 5.3 3.5 2.1 1.5 1 0.7]'; % amplitude of each SS function in order in PControl
name = 'amp'; % name of identifier at end of file name

%% EXPERIMENTAL PARAMETERS
n_tracktime = 10 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 1;                  % pattern ID
yPos = 5;                   % 30 deg spatial frequency
xUpdate = 400;           	% function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 0:2;                	% Analog input channels
AO = 1;                     % Analog output channels

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
% [s,~] = MC_USB_1208FS_PLUS(Fs,AI,AO);
[s,~] = NI_USB_6212(Fs,AI,AO);

%% Camera Trigger Signal
off = 0.1;
t = 0:1/s.Rate:n_tracktime + off;
TRIG = ((1/2)*(square(2*pi*FPS*t,5) - 1)');
TRIG(TRIG==-1) = 4;
end_off = round(Fs*off);
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
% [vid,src] = Basler_acA640_750um(nFrame);
[vid,~] = Basler_acA640_120gm(FPS,Gain,nFrame);

%% Set variable to control positionn function
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
freq_all = freq(func(func_all)); % true values (frequency)
n_trial = n_rep * n_func;

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   Amp = %i \n', ii, val_all(ii))
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
    Panel_com('send_gain_bias', [0 0 0 0])
	
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
    
    % CL BUFFER
    Arena_CL(1,'X',-15)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_' name '_' ...
        num2str(val_all(ii)) '_freq_' num2str(freq_all(ii)) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end