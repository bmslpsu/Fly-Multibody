function [] = Experiment_SS_amp(Fn)
%% Experiment_SOS: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 100;
%% Set directories & experimental parameters
root = 'C:\BC\Experiment_SS_norm_250';

%% EXPERIMENTAL PARAMETERS
n_tracktime = 11;           % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 10;                 % # of repetitions
patID = 1;                  % pattern ID
yPos = 5;                   % 30 deg spatial frequency
xUpdate = 500;              % function update rate
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

%% Set variable to control pattern oscillation frequency
amp = [60 26.25 3.75 41.25 18.75 11.25 7.5];  % make sure same order as on SD card
freq = [0.7 1.5 10.6 1 2.1 3.5 5.3]'; % make sure same order as on SD card
n_freq = length(freq); % # of frequencies
freqI = (1:n_freq)'; % oscillation frequency indicies

% Create sequence of randomly shuffled frequencies
freqI_all = nan(n_freq*n_rep,1);
pp = 0;
for kk = 1:n_rep
    freq_rand = freqI(randperm(n_freq),:);    % reshuffle randomly
    freqI_all(pp+1:pp+n_freq,1) = freq_rand;  % add rep
    pp = kk*n_freq;
end

freq_all = freqI_all;
amp_all = freqI_all;
for kk = 1:length(freq)
    freq_all(freqI_all==freqI(kk)) = freq(kk);
    amp_all(freqI_all==freqI(kk)) = amp(kk);
end

n_trial = n_rep * n_freq;

%% EXPERIMENT LOOP
disp('Start Experiment:')

for ii = 1:n_trial
    fprintf('Trial: %i  Freq: %1.2f Amp: %1.2f \n', ii, freq_all(ii), amp_all(ii))
    preview(vid) % open video preview window
    
    Panel_com('stop')
    pause(2) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    Panel_com('set_pattern_id', patID);	% set pattern
    pause(n_pause)
    Panel_com('set_position', [1, yPos]); % set starting position (xpos,ypos)
    pause(n_pause)
    Panel_com('set_posfunc_id',[1,freqI_all(ii)]); % arg1 = channel (x=1,y=2); arg2 = funcID
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
    Arena_Ramp(1,16)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_freq_' ...
        num2str(freq_all(ii)) '_Amp_' num2str(amp_all(ii)) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
end

delete(vid)
disp('Done');
daqreset
imaqreset
end