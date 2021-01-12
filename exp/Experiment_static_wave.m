function [] = Experiment_static_wave(Fn)
%% Experiment_static_wave: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_static_wave';
name = 'wave'; % name of identifier at end of file name

%% EXPERIMENTAL PARAMETERS
n_tracktime = 20 + 1;     	% length(func)/fps seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 1;                  % pattern ID
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
[vid,src] = Basler_acA640_750um(nFrame);

%% Set variable to control y-position function
val = 7.5*[0 3 4 8 inf]';   % [°] wavelengths
n_val = length(val);       	% # of values
ypos = [1 4 5 7 12];    	% pattern y-pos corresponding to each wavelength

% Create sequence of randomly shuffled frequencies
val_all = nan(n_val*n_rep,1);
pp = 0;
for kk = 1:n_rep
    val_rand = val(randperm(n_val),:);    % reshuffle randomly
    val_all(pp+1:pp+n_val,1) = val_rand;  % add rep
    pp = kk*n_val;
end

% y-pos index vector
ypos_all = val_all;
for kk = 1:length(ypos)
    ypos_all(ypos_all==val(kk)) = ypos(kk);
end

n_trial = n_rep * n_val;

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 1:n_trial
    fprintf('Trial: %i   val = %f \n', ii, val_all(ii))
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
    Panel_com('set_position', [randi(96), ypos_all(ii)]); % set starting position (xpos,ypos)
    pause(n_pause)
    Panel_com('set_posfunc_id',[1,1]); % arg1 = channel (x=1,y=2); arg2 = funcID
    pause(n_pause)
	Panel_com('set_funcX_freq', 50); % update rate for x-channel
    pause(n_pause)
    Panel_com('set_funcY_freq', 50); % update rate for y-channel
    pause(n_pause)
    Panel_com('set_mode', [0,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
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
    
    % SPIN BUFFER
    Arena_Ramp(1,10)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_' name '_' ...
        num2str(val_all(ii)) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end