function [] = Spin(Fn,root,gain)
%% Spin: runs a spin trial for the magnetic tether
% Fn is the fly number
daqreset
imaqreset

spinroot = fullfile(root,'spin');
vel = 3.75*gain; % spin velocity [deg/s]

% PARAMETERS
n_tracktime = 11;       	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_trial = 1;                % # of repetitions
patID = 2;                  % Spatial frequency grating pattern
yPos  = 7;                  % 30 deg spatial frequency
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Gain = 11.993167134727036;	% camera gain
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 0:2;                	% Analog input channels
AO = 1;                     % Analog output channels

%% Set up data acquisition on NiDAQ (session mode)
% DAQ Setup
[s,~] = MC_USB_1208FS_PLUS(Fs,AI,AO);

% Camera Trigger Signal
t = 0:1/s.Rate:n_tracktime;
TRIG = ((1/2)*(square(2*pi*FPS*t,5) - 1)');
TRIG(TRIG==-1) = 4;

% Camera Setup
[vid,src] = Basler_acA640_750um(FPS,Gain,nFrame);

%% SPIN LOOP
disp('Start Spin:')
for ii = 1:n_trial   
    disp('Spin Trial')
    disp(num2str(ii));  % print counter to command line
    preview(vid);       % open video preview window
	    
    pause(1)
    
    % EXPERIMENT SETUP
    disp('Play Stimulus: ')
    Panel_com('set_pattern_id', patID); pause(n_pause)          % set pattern
    Panel_com('set_position',   [15 , yPos]); pause(n_pause) 	% set starting position (xpos,ypos)
	Panel_com('set_funcX_freq', 50); pause(n_pause)             % update rate for x-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)             % update rate for y-channel
    Panel_com('set_mode',       [0,0]); pause(n_pause)         	% 0=open,1=closed,2=fgen,3=vmode,4=pmode
	Panel_com('send_gain_bias',	[gain,0,0,0]); pause(n_pause) 	% set gain

    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    queueOutputData(s,TRIG) % set trigger AO signal
    T = timer('StartDelay',0.5,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        Panel_com('stop') % stop stimulus
        stop(vid) % stop video buffer
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    pause(1)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_vel_' num2str(vel) '_Spin.mat'];
    save(fullfile(spinroot,fname),'-v7.3','data','t_p','vidData','t_v');
end

delete(vid)
disp('Done')
daqreset
imaqreset
end