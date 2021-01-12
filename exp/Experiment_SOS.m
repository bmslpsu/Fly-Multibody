function [] = Experiment_SOS(Fn)
%% Experiment_SOS: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_SOS_v2';
% Fn = 200;

% Spin Trial
% Pick random direction
% dir = 0;
% while dir==0
%     dir = randi([-1,1],1);
% end
% Spin(Fn,root,dir*16);

%% EXPERIMENTAL PARAMETERS
n_tracktime = 21;           % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_trial = 20;               % # of repetitions
patID = 1;                  % Spatial frequency grating pattern
yPos  = 5;                  % 30 deg spatial frequency
funcX = 1;                  % SOS (20s)
xUpdate = 200;              % function update rate
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
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
[vid,src] = Basler_acA640_750um(nFrame);

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 1:n_trial   
    disp('Trial')
    disp(num2str(ii));  % print counter to command line
    preview(vid);       % open video preview window
	
    % SPIN BUFFER
    Arena_Ramp(2,16)
    pause(n_resttime)
    Panel_com('stop')
    
    pause(1) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    Panel_com('set_pattern_id', patID);	% set pattern
    pause(n_pause)
    Panel_com('set_position', [randi(96), yPos]); % set starting position (xpos,ypos)
    pause(n_pause)
    Panel_com('set_posfunc_id',[funcX, 1]); % arg1 = channel (x=1,y=2); arg2 = funcID
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
    Arena_Ramp(2,16)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_SOS.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
end

delete(vid)
disp('Done');
daqreset
imaqreset
end