function [] = make_SOS_roll()
%% make_SOS: make the dataset for head-free/body-free flies in the magnetic tether for SOS stimuli
%

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2_roll';
savedir = 'E:\DATA\Multisensory_control';

% Set state variable to position or velocity
clss = 'position';
% clss = 'velocity';

% Get experiment name
[~,exp_name,~] = fileparts(rootdir);
exp_name = textscan(char(exp_name), '%s', 'delimiter', '_');
exp_name = string(exp_name{1});
exp_name = char(strjoin(exp_name(2:end), '_'));
filename = [exp_name '_' num2str(clss)];

%% Setup Directories
root.base = rootdir;
root.reg = fullfile(root.base,'registered');
root.body = fullfile(root.reg,'reg_angles');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');
root.head_roll = fullfile(root.head ,'tracked_head_roll');
root.func = fullfile(root.base ,'function');

% Load function files
func_list = dir(root.func);
func_list = func_list(~[func_list.isdir]);
n_cond = length(func_list);
FUNC = cell(n_cond,1);
for f = 1:n_cond
    FUNC{f} = load(fullfile(root.func, func_list(f).name));
    FUNC{f}.name = func_list(f).name;
end

% Select files
[D,I,N,U,T,~,~,basename] = get_file_data(root.head_roll,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = get_file_data(root.body,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = get_file_data(root.benifly,'*.csv',false);

%% Get Data
close all
clc
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Body saccade detection parameters
scd.threshold = [20, 1, 3, 0];
scd.true_thresh = 220;
scd.Fc_detect = [40 nan];
scd.Fc_ss = [20 nan];
scd.amp_cut = 5;
scd.dur_cut = 1;
scd.direction = 0;
scd.direction = 0;
scd.pks = [];
scd.sacd_length = nan;
scd.min_pkdist = 0.2;
scd.min_pkwidth = 0.02;
scd.min_pkprom = 50;
scd.min_pkthresh = 0;
scd.boundThresh = [0.2 40];

Fs = 100;
Fc = 40;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
showplot = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N{1,3});
DATA = [D , splitvars(table(num2cell(zeros(N.file,6))))];
DATA.Properties.VariableNames(4:end) = {'reference', 'body', 'head_yaw', 'head_roll', 'error', 'body_saccade'};
for n = 1:N.file
    disp(basename{n})
    
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'angles'); % load body angles
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
    data.roll = load(fullfile(root.head_roll, [basename{n} '.mat']),'roll'); % load head roll angles
    
    % Get synced frame times and pattern data
    daq_time = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger = data.daq.data(:,1);
    [TRIG,PAT] = sync_pattern_trigger(daq_time, daq_pattern, trigger, func_length, [],  false);
    trig_time = TRIG.time_sync;
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.angles;
	head_yaw = hampel(trig_time, data.head.head_data.angle, 3);
    head_roll = hampel(trig_time, data.roll.roll, 20);
    
    % Interpolate so all signals have the same times
    %[b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    %Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    Reference = 3.75*(round(Reference/3.75));
    
    Body = interp1(trig_time, body,  tintrp, 'pchip');
    Head_yaw = interp1(trig_time, head_yaw,  tintrp, 'pchip');
    Head_roll = interp1(trig_time, head_roll,  tintrp, 'pchip');
    
    % Detect & remove saccades
    body_scd = saccade(Body, tintrp, scd, false);
    
    % Store signals
    n_detrend = 5;
    DATA.body_saccade{n}    = body_scd;
    DATA.reference{n}       = singal_attributes(Reference, tintrp);
    DATA.body{n}            = singal_attributes(body_scd.shift.IntrpPosition, tintrp, Fc, n_detrend);
    DATA.head_yaw{n}       	= singal_attributes(Head_yaw, tintrp, Fc, 1);
    DATA.head_roll{n}     	= singal_attributes(Head_roll, tintrp, Fc, 1);
    
    Error                   = DATA.reference{n}.position - DATA.body{n}.position - DATA.head_yaw{n}.position;
	DATA.error{n}           = singal_attributes(Error, tintrp, Fc);

    % Debug plot
    if showplot
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            plot(tintrp, DATA.head_yaw{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, 5*DATA.dwba{n}.position, 'm', 'LineWidth', 1)
            leg = legend('Reference','Body','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
           	
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = sort(FUNC{I{n,3}}.All.Freq, 'ascend');
    REF = DATA.reference{n}.(clss);
    BODY = DATA.body{n}.(clss);
    HEAD_YAW = DATA.head_yaw{n}.(clss);
    HEAD_ROLL = DATA.head_roll{n}.(clss);
    ERROR = DATA.error{n}.(clss);
    
    SYS_ref2_head_body = frf(tintrp, REF , IOFreq, false, BODY, HEAD_YAW);
	SYS_body2_head = frf(tintrp, BODY, IOFreq, false, HEAD_YAW);
 	SYS_err2_head_body = frf(tintrp, ERROR, IOFreq, false, BODY, HEAD_YAW);
    SYS_yaw2roll = frf(tintrp, HEAD_YAW, IOFreq, false, HEAD_ROLL);
    
	SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_body2_head, SYS_err2_head_body, SYS_yaw2roll);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_all;
end

%% Group Data
[GRAND, ~] = grand_stats(ALL);

%% Plot
yaw = cellfun(@(x) x.position, DATA.head_yaw, 'UniformOutput', false);
yaw = cat(2, yaw{:});
roll = cellfun(@(x) x.position, DATA.head_roll, 'UniformOutput', false);
roll = cat(2, roll{:});

yaw_color = [0 0.2 0.8];
roll_color = [0.8 0.1 0.8];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 10 2.5])

clear ax
ax(1) = subplot(1,5,1:3); cla ; hold on
     plot(tintrp, mean(yaw,2), 'Color', [yaw_color 1], 'LineWidth', 0.75)
     plot(tintrp, mean(roll,2), 'Color', [roll_color 1], 'LineWidth', 0.75)
     xlim([-0.5 20])
     ylim(5*[-1 1])
     xlabel('Time (s)')
     ylabel('Head displacement (°)')
     
ax(2) = subplot(1,5,4); cla ; hold on
    bin_size = 0.5;
    bin_range = 10;
    bins = -(bin_range + bin_size/2):bin_size:(bin_range + bin_size/2);
    histogram(yaw(:), bins, 'FaceColor', yaw_color, 'EdgeColor', 'none', 'Normalization', 'probability')
    histogram(roll(:), bins, 'FaceColor', roll_color, 'EdgeColor', 'none', 'Normalization', 'probability')
  	xlim(bin_range*[-1 1])
    ylabel('Probability')
    xlabel('Head displacement (°)')
    
ax(3) = subplot(1,5,5); cla ; hold on
     %plot(GRAND.fly_stats.mean.Fv.mean(:,1), GRAND.fly_stats.mean.Cohr.mean(:,4), ...
         %'Color', 'k', 'LineWidth', 0.75)
     [~,h.line] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,4), GRAND.fly_stats.mean.IOCohr.std(:,4), ...
         GRAND.fly_stats.mean.IOFv.mean(:,1), 1, 1, 'k', 'k', 0.2, 0.75);
     set(h.line, 'Marker', '.', 'MarkerSize', 12, 'MarkerFaceColor', 'none')
     xlim([0.2 15])
     xticks([0.2 1 10])
     ylim([0 1])
     set(ax(3), 'XScale', 'log')
     xlabel('Frequency (hz)')
     ylabel('Coherence')
     
set(ax, 'Color', 'none', 'LineWidth', 0.75)

%% SAVE
disp('Saving...')
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')

end