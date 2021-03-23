function [] = make_head_free_magno_static(rootdir)
%% make_head_free_magno_static:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_static_wave';
% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_static_wave_head_fixed';

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.body = fullfile(root.base,'tracked_body');
root.reg = fullfile(root.base,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');
head_test = exist(root.head, 'file');

if head_test == 7
    filename = 'Static_HeadFree';
    
    % Select files
    [D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
    % [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);
else
	filename = 'Static_HeadFixed';
    
	% Select files
    [D,I,N,U,T,~,~,basename] = GetFileData(root.body,'*.mat',false);
end

%% Get Data
close all
clc

% Body saccade detection parameters
body_scd.threshold = [20, 1, 3, 0];
body_scd.true_thresh = 220;
body_scd.Fc_detect = [40 nan];
body_scd.Fc_ss = [20 nan];
body_scd.amp_cut = 5;
body_scd.dur_cut = 1;
body_scd.direction = 0;
body_scd.direction = 0;
body_scd.pks = [];
body_scd.sacd_length = nan;
body_scd.min_pkdist = 0.2;
body_scd.min_pkwidth = 0.02;
body_scd.min_pkprom = 50;
body_scd.min_pkthresh = 0;
body_scd.boundThresh = [0.2 40];

head_scd.threshold = [40, 1, 2, 0];
head_scd.true_thresh = 200;
head_scd.Fc_detect = [15 nan];
head_scd.Fc_ss = [nan nan];
head_scd.amp_cut = 4;
head_scd.dur_cut = 0.1;
head_scd.direction = 0;
head_scd.direction = 0;
head_scd.pks = [];
head_scd.sacd_length = nan;
head_scd.min_pkdist = 0.2;
head_scd.min_pkwidth = 0.02;
head_scd.min_pkprom = 50;
head_scd.min_pkthresh = 0;
head_scd.boundThresh = [0.2 60];

Fs = 100;
Fc = 40;
func_length = 20;
startI = round(5000*0.5);
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
DATA = [D , splitvars(table(num2cell(zeros(N.file,11))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','body_raw','head','error',...
    'dwba','lwing','rwing','body2head','body_saccade','head_saccade'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...  % load head & wing angles from Benifly
                            [basename{n} '.csv']));
    if head_test == 7
        data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
    else % no head data
        data.head.head_data.angle = 0*data.body.bAngles;
    end
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,~]    = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, startI, false, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    lwing = rad2deg(data.benifly.LWing);
    rwing = rad2deg(data.benifly.RWing);
    lwing = hampel(data.benifly.Time, lwing);
    rwing = hampel(data.benifly.Time, rwing);
	lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
    % Get head & body anlges
    body = data.body.bAngles;
	head = data.head.head_data.angle;
    
    % Interpolate so all signals have the same times    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Detect & remove saccades
    body_saccade = saccade_v1(Body, tintrp, body_scd, false);
    head_saccade = saccade_v1(Head, tintrp, head_scd, false);
    
    %figure (1) ; pause ; close all

    % Store signals
    n_detrend = 5;
    DATA.body_saccade{n}    = body_saccade;
    DATA.head_saccade{n}    = head_saccade;
    DATA.reference{n}       = nan;
    DATA.body{n}            = singal_attributes(body_saccade.shift.IntrpPosition, tintrp, [], n_detrend);
    DATA.body_raw{n}      	= singal_attributes(Body, tintrp);
    DATA.head{n}            = singal_attributes(Head, tintrp, []);
    DATA.error{n}           = nan;
    DATA.dwba{n}            = singal_attributes(dWBA, tintrp, [], n_detrend);
    DATA.lwing{n}           = singal_attributes(LWing, tintrp, [], n_detrend);
    DATA.rwing{n}           = singal_attributes(RWing, tintrp, [], n_detrend);
    
    body2head.coherence = mscohere(DATA.body{n}.position, DATA.head{n}.position, ...
                                            [], [], DATA.body{n}.Fv, Fs);
    DATA.body2head{n} = body2head;
    
    
    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            %plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            %plot(tintrp, DATA.dwba{n}.position, 'm', 'LineWidth', 1)
            %leg = legend('Body','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
          	
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        %linkaxes(ax,'x')
        pause
    end
end

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end