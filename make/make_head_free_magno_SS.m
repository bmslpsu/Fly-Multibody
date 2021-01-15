function [] = make_head_free_magno_SS(rootdir)
%% make_head_free_magno_SS:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250';
exp_name = textscan(char(rootdir), '%s', 'delimiter', '_');
exp_typ = exp_name{1}{end-1}; % type of stimuli (vel or pos)
exp_ver = exp_name{1}{end}; % version of experiment (v1, v2, ...)

clss = 'position';
% clss = 'velocity';
filename = ['SS_HeadFree_' exp_typ '_' exp_ver '_' num2str(clss)];

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.body = fullfile(root.base,'tracked_body');
root.reg = fullfile(root.base,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');
root.func = fullfile(root.base ,'function');

% Load function files
func_list = dir(root.func);
func_list = func_list(~[func_list.isdir]);
n_cond = length(func_list); % number of stimuli (functions) used in experiment
FUNC = cell(n_cond,1);
for f = 1:n_cond
    FUNC{f} = load(fullfile(root.func, func_list(f).name));
    FUNC{f}.name = func_list(f).name;
end
forder = [4 7 1 2 3 5 6];
FUNC = FUNC(forder,1);

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

%% Get Data %%
close all
clc

scd.thresh = [20, 1, 3, 0];
scd.true_thresh = 220;
scd.Fc_detect = [40 nan];
scd.Fc_ss = [20 nan];
scd.amp_cut = 5;
scd.dur_cut = 1;
scd.direction = 0;
scd.direction = 0;
scd.pks = [];
scd.sacd_length = nan;
scd.min_pkdist = 0.5;
scd.min_pkwidth = 0.02;
scd.min_pkprom = 50;
scd.min_pkthresh = 0;
scd.boundThresh = [0.2 40];

Fs = 100;
Fc = 40;
func_length = 10;
tintrp = (0:(1/Fs):func_length)';
debug = true;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N{1,3});
DATA = [I , splitvars(table(num2cell(zeros(N.file,7))))];
DATA.Properties.VariableNames(5:end) = {'reference','body','head','error','dwba','lwing','rwing'};
for n = 53:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
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
    
  	% Filter wing angles
    lwing = rad2deg(data.benifly.LWing);
    rwing = rad2deg(data.benifly.RWing);
    lwing = hampel(data.benifly.Time, lwing);
    rwing = hampel(data.benifly.Time, rwing);
	lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.bAngles;
	head = data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    Fc_pat_ratio = 2;
    %[b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    [b_pat, a_pat] = butter(3, Fc_pat_ratio * D.freq(n) / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = Reference - mean(Reference);
    Reference = 3.75*round(Reference/3.75);
    %Reference = filtfilt(b_pat, a_pat, Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Body - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Detect & remove saccades
    body_scd = saccade_v1(Body, tintrp, scd.thresh, scd.true_thresh, scd.Fc_detect, ...
                            scd.Fc_ss, scd.amp_cut, scd.dur_cut , scd.direction, scd.pks, ...
                            scd.sacd_length, scd.min_pkdist, scd.min_pkwidth, scd.min_pkprom, ...
                            scd.min_pkthresh, scd.boundThresh, false);
%     figure (1)
%     pause
%     close all
    
    % Store signals
    n_detrend = 7;
    DATA.reference{n}   = singal_attributes(Reference, tintrp);
    DATA.body{n}        = singal_attributes(body_scd.shift.IntrpPosition, tintrp, [], n_detrend);
    DATA.body{n}        = singal_attributes(DATA.body{n}.detrend, tintrp, [], []);
    DATA.head{n}        = singal_attributes(Head, tintrp, [], n_detrend);
    DATA.head{n}        = singal_attributes(DATA.head{n}.detrend, tintrp, [], []);
    DATA.error{n}       = singal_attributes(Error, tintrp, [], n_detrend);
    DATA.dwba{n}    	= singal_attributes(dWBA, tintrp, 1.5 * D.freq(n), n_detrend);
    DATA.dwba{n}        = singal_attributes(DATA.dwba{n}.detrend, tintrp, [], []);
    DATA.dwba{n}        = singal_attributes(DATA.dwba{n}.position, tintrp, 1.5 * D.freq(n), []);
    
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, []);
    DATA.rwing{n}       = singal_attributes(RWing, tintrp, []);
    
%     hold on
%     plot(tintrp, DATA.body{n}.trend, 'g--', 'LineWidth', 1)
%     plot(tintrp, Body, 'k', 'LineWidth', 1)
%     plot(tintrp, DATA.body{n}.position, 'b', 'LineWidth', 1)
%     plot(tintrp, DATA.body{n}.detrend, 'r', 'LineWidth', 1)
%     pause
%     cla

    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, 5*DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Reference','Body','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
%     hold on
%     plot(DATA.reference{n}.time, DATA.reference{n}.position)
%     plot(FUNC{I{n,3}}.All.time, FUNC{I{n,3}}.All.X, 'r')
%     plot(FUNC{I{n,3}}.All.time, FUNC{I{n,3}}.All.X_step, 'g')
%     pause
%     cla
    
    IOFreq = sort(FUNC{I{n,3}}.All.Freq, 'ascend');
    REF = DATA.reference{n}.(clss);
    BODY = DATA.body{n}.(clss);
    HEAD = DATA.head{n}.(clss);
    dWBA = DATA.dwba{n}.(clss);
    LWING = DATA.lwing{n}.(clss);
    RWING = DATA.rwing{n}.(clss);
    
    SYS_ref2_head_body = frf(tintrp, REF , IOFreq, false, BODY, HEAD);
    SYS_ref2_wing = frf(tintrp, REF, IOFreq, false, LWING, RWING, dWBA);
    SYS_head2_body_wing = frf(tintrp, HEAD, IOFreq, false, BODY, dWBA);
    SYS_wing2_body = frf(tintrp, dWBA, IOFreq, false, BODY);
    %SYS_left2_right = frf(tintrp, LWING, IOFreq, false, RWING);
    
    SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_ref2_wing, ...
                                    SYS_head2_body_wing, SYS_wing2_body);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_all;
end

%% Group Data
% clc
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N{1,3}
    GRAND.all(v) = cell2struct(cell(nfield,1),fields);
    GRAND.all_trial(v) = cell2struct(cell(nfield,1),fields);
    %GRAND.fly_all(fr) = cell2struct(cell(n_stat_fields,1),stat_fields);
    %GRAND.fly_all(fr) = [];
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
ax = gobjects(N.amp,1);
cc = hsv(N.amp);
pp = 1;
for a = fliplr(1:N.amp)
    ax(a) = subplot(N.amp,1,pp); cla ; hold on
    plot(FUNC{a}.All.time, FUNC{a}.All.X, 'LineWidth', 1, 'Color', 'k')
        %plot(squeeze(GRAND.all(a).Time), squeeze(GRAND.all(a).State(:,1,:)), ...
            %'LineWidth', 0.25, 'Color', [0.5 0.5 0.5 0.3])
%     plot(median(squeeze(GRAND.all(a).Time),2), median(squeeze(GRAND.all(a).State(:,1,:)),2), ...
%         'LineWidth', 1, 'Color', [cc(a,:) 1])
    [h1,h2] = PlotPatch(10*median(squeeze(GRAND.all(a).State(:,4,:)),2), ...
        10*std(squeeze(GRAND.all(a).State(:,4,:)),[],2), ...
        median(squeeze(GRAND.all(a).Time),2),...
        1,1,cc(a,:), 0.7*cc(a,:), 0.2, 2);
    pp = pp + 1;
end

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end