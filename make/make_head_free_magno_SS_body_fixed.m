function [] = make_head_free_magno_SS_body_fixed(rootdir)
%% make_head_free_magno_SS_body_fixed:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_amp_3.75_body_fixed';
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250_body_fixed';
exp_name = textscan(char(rootdir), '%s', 'delimiter', '_');
exp_typ = exp_name{1}{end-3}; % type of stimuli (vel or pos)
exp_ver = exp_name{1}{end-2}; % version of experiment (v1, v2, ...)

clss = 'position';
% clss = 'velocity';
filename = ['SS_BodyFixed_' exp_typ '_' exp_ver '_' num2str(clss)];

%% Setup Directories %%
root.base = rootdir;
root.benifly = fullfile(root.base ,'tracked_head_wing');
root.head = fullfile(root.base ,'tracked_head_tip');
root.func = fullfile(root.base ,'function');

% Load function files
func_list = dir(root.func);
func_list = func_list(~[func_list.isdir]);
n_cond = length(func_list); % number of stimuli (functions) used in experiment
FUNC = cell(n_cond,1);
freq_order = nan(n_cond,1);
for f = 1:n_cond
    FUNC{f} = load(fullfile(root.func, func_list(f).name));
    FUNC{f}.name = func_list(f).name;
    freqI = strfind(func_list(f).name, 'freq');
    matI = strfind(func_list(f).name, 'mat');
    freq_order(f) = str2double(func_list(f).name(freqI+5:matI-2));
end
[~,forder] = sort(freq_order);
FUNC = FUNC(forder,1);

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

%% Get Data %%
close all
clc

Fs = 100;
Fc = 40;
func_length = 10;
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N{1,4});
DATA = [D , splitvars(table(num2cell(zeros(N.file,8))))];
DATA.Properties.VariableNames(5:end) = {'reference','body','head','error',...
    'dwba','lwing','rwing','body_saccade'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...  % load head & wing angles from Benifly
                            [basename{n} '.csv']));
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, [], false, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles (flipped)
    lwing = rad2deg(data.benifly.RWing);
    rwing = rad2deg(data.benifly.LWing);
    lwing = hampel(data.benifly.Time, lwing);
    rwing = hampel(data.benifly.Time, rwing);
	lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
    % Get pattern, head, & body anlges
    pat = 3.75*unwrap_custom(PAT.pos, 50, 1, 96);
    %pat = 3.75*PAT.pos;
	head = -data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    Fc_pat_ratio = 2;
    %[b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    [b_pat, a_pat] = butter(3, Fc_pat_ratio * D.freq(n) / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = Reference - mean(Reference);
    Reference = 3.75*round(Reference/3.75);
    %Reference = filtfilt(b_pat, a_pat, Reference);
    
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    fc_wing = 1.5 * D.freq(n);
    n_detrend = 5;
    DATA.reference{n}	= singal_attributes(Reference, tintrp);
    DATA.head{n}       	= singal_attributes(Head, tintrp, [], n_detrend);
    DATA.error{n}     	= singal_attributes(Error, tintrp, [], n_detrend);
    DATA.dwba{n}       	= singal_attributes(dWBA, tintrp, fc_wing, n_detrend);
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, fc_wing, n_detrend);
    DATA.rwing{n}    	= singal_attributes(RWing, tintrp, fc_wing, n_detrend);

    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 0.5)
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, DATA.dwba{n}.position, 'm', 'LineWidth', 1)
            leg = legend('Reference','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = sort(FUNC{I.freq(n)}.All.Freq, 'ascend');
    REF = DATA.reference{n}.(clss);
    HEAD = DATA.head{n}.(clss);
    ERROR = DATA.error{n}.(clss);
    dWBA = DATA.dwba{n}.position;
    LWING = DATA.lwing{n}.position;
    RWING = DATA.rwing{n}.position;
    
    SYS_ref2_head = frf(tintrp, REF , IOFreq, false, HEAD);
    SYS_ref2_wing = frf(tintrp, REF, IOFreq, false, dWBA);
    SYS_head2wing = frf(tintrp, HEAD, IOFreq, false, dWBA);
    SYS_ref2_left_right = frf(tintrp, REF, IOFreq, false, LWING, RWING);
    
    SYS_all = CatStructFields(2, SYS_ref2_head, SYS_ref2_wing, SYS_head2wing, SYS_ref2_left_right);
    
    ALL{I.fly(n),I.freq(n)}(end+1,1) = SYS_all;
end

%% Group Data
% clc
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N.freq
    GRAND.all(v) = cell2struct(cell(nfield,1),fields);
    GRAND.all_trial(v) = cell2struct(cell(nfield,1),fields);
    for n = 1:N.fly
        for f = 1:nfield
            FLY.all(n,v).(fields{f})    = cat(3,ALL{n,v}.(fields{f}));
            FLY.stats(n,v).(fields{f})  = system_stats(FLY.all(n,v).(fields{f}),3);
            GRAND.all(v).(fields{f}) 	= cat(3,GRAND.all(v).(fields{f}), ...
                                                        FLY.all(n,v).(fields{f}));

            stat_fields = fieldnames(FLY.stats(n,v).(fields{f}));
            n_stat_fields = length(stat_fields);
            for s = 1:n_stat_fields
                GRAND.fly_all(v).(stat_fields{s}).(fields{f})(:,:,n) = ...
                    FLY.stats(n,v).(fields{f}).(stat_fields{s});
            end
        end
    end
    
    for f = 1:nfield
        stat_fields = fieldnames(FLY.stats(n,v).(fields{f}));
        n_stat_fields = length(stat_fields);
        for s = 1:n_stat_fields
            GRAND.fly_stats(v).(stat_fields{s}).(fields{f}) = ...
                system_stats(GRAND.fly_all(v).(stat_fields{s}).(fields{f}),3);
        end
    end
    GRAND.all_trial(v) = structfun(@(x) system_stats(x,3), GRAND.all(v), 'UniformOutput', false);
end

%%
fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches')
ax = gobjects(N{1,3},1);
cc = hsv(N{1,3});
pp = 1;
for a = 1:N.freq
    ax(a) = subplot(N.freq,1,pp); cla ; hold on
    plot(FUNC{a}.All.time, FUNC{a}.All.X, 'LineWidth', 1, 'Color', 'k')
        %plot(squeeze(GRAND.all(a).Time), squeeze(GRAND.all(a).State(:,1,:)), ...
            %'LineWidth', 0.25, 'Color', [0.5 0.5 0.5 0.3])
%     plot(median(squeeze(GRAND.all(a).Time),2), median(squeeze(GRAND.all(a).State(:,1,:)),2), ...
%         'LineWidth', 1, 'Color', [cc(a,:) 1])
    scl = median(abs(median(squeeze(GRAND.all(a).State(:,1,:)),2))) / median(abs(median(squeeze(GRAND.all(a).State(:,2,:)),2)));
    [h1,h2] = PlotPatch(median(squeeze(GRAND.all(a).State(:,1,:)),2), ...
        std(squeeze(GRAND.all(a).State(:,5,:)),[],2), ...
        median(squeeze(GRAND.all(a).Time),2),...
        0, 1, [0 0 1], 0.7*[0 0 1], 0.2, 1);
    [h1,h2] = PlotPatch(median(squeeze(GRAND.all(a).State(:,2,:)),2), ...
        std(squeeze(GRAND.all(a).State(:,1,:)),[],2), ...
        median(squeeze(GRAND.all(a).Time),2),...
        0, 1, [1 0 0], 0.7*[1 0 0], 0.2, 1);
    pp = pp + 1;
end
% linkaxes(ax)

%%
fig = figure (2); clf
set(fig, 'Color', 'w', 'Units', 'inches')
ax = gobjects(3,1);
cc = hsv(N.freq);
pp = 1;
for a = fliplr(1:N.freq)
    ax(1) = subplot(3,1,1) ; hold on
        plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOGain.mean(1), '*b')
      	plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOGain.mean(2), '*r')
        plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOGain.mean(3), '*m')

    ax(2) = subplot(3,1,2) ; hold on
        plot(GRAND.fly_stats(a).mean.IOFv.mean, rad2deg(GRAND.fly_stats(a).circ_mean.IOPhaseDiff.circ_mean(1)), '*b')
      	plot(GRAND.fly_stats(a).mean.IOFv.mean, rad2deg(GRAND.fly_stats(a).circ_mean.IOPhaseDiff.circ_mean(2)), '*r')
        plot(GRAND.fly_stats(a).mean.IOFv.mean, rad2deg(GRAND.fly_stats(a).circ_mean.IOPhaseDiff.circ_mean(3)), '*m')
        
    ax(3) = subplot(3,1,3) ; hold on
    	plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOCohr.mean(:,1), '*b')
      	plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOCohr.mean(:,2), '*r')
        plot(GRAND.fly_stats(a).mean.IOFv.mean, GRAND.fly_stats(a).mean.IOCohr.mean(:,3), '*m')
        
    pp = pp + 1;
end
set(ax, 'LineWidth', 1.5)
set(ax(1), 'YLim', [0 1])
set(ax(3), 'YLim', [0 1])
% linkaxes(ax)

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end