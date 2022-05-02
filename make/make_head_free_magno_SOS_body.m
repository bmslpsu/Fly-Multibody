function [] = make_head_free_magno_SOS_body()
%% make_head_free_magno_SOS_body:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_52_add_mass';
set_mass = 3300;
root_set = fullfile(rootdir, num2str(set_mass));

% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_52';
% set_mass = 100;
% root_set = rootdir;

warning('off', 'signal:findpeaks:largeMinPeakHeight')

clss = 'position';
clss = 'velocity';

filename = ['SOS_' num2str(set_mass) '_' num2str(clss)];

%% Setup Directories
root.base = root_set;
root.body = fullfile(root.base,'tracked_body');
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
[D,I,N,U,T,~,~,basename] = GetFileData(root.body,'*.mat',false);

%% Get Data
close all
clc

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
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N{1,3});
DATA = [D , splitvars(table(num2cell(zeros(N.file,8))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error',...
    'dwba','lwing','rwing','body_saccade'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	%data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles', 'data', 't_p'); % load body angles
    
    % Get synced frame times and pattern data
    daq_time    = data.body.t_p;
    daq_pattern = data.body.data(:,2);
    trigger     = data.body.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, [], false, false);
    trig_time   = TRIG.time_sync;
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.bAngles;
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    %Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    Reference = 3.75*(round(Reference/3.75));
    
    Body = interp1(trig_time, body,  tintrp, 'pchip');
    Body = Body - mean(Body);
    
    % Detect & remove saccades
    body_scd = saccade_v1(Body, tintrp, scd, false);
    
    % Store signals
    n_detrend = 5;
    DATA.body_saccade{n}    = body_scd;
    DATA.reference{n}       = singal_attributes(Reference, tintrp);
    DATA.body{n}            = singal_attributes(body_scd.shift.IntrpPosition, tintrp, [], n_detrend);
    
    Error                   = DATA.reference{n}.position - DATA.body{n}.position;
	DATA.error{n}           = singal_attributes(Error, tintrp);

%     hold on
%     plot(tintrp, body_scd.shift.IntrpPosition, 'b', 'LineWidth', 1)
%     plot(tintrp, DATA.body{n}.trend, 'g--', 'LineWidth', 1)
%     plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
%     %plot(tintrp, Head, 'k', 'LineWidth', 1)
%     %plot(tintrp, DATA.head{n}.trend, 'g--', 'LineWidth', 1)
%     %plot(tintrp, DATA.head{n}.position, 'r', 'LineWidth', 1)
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
           	
       	set_mass(gcf, 'Color', 'w')
       	set_mass(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = sort(FUNC{I{n,3}}.All.Freq, 'ascend');
    IOFreq = IOFreq(1:end-1);
    REF = DATA.reference{n}.(clss);
    BODY = DATA.body{n}.(clss);
    ERROR = DATA.error{n}.(clss);

    SYS_ref2_body = frf(tintrp, REF , IOFreq, false, BODY);
 	SYS_err2_body = frf(tintrp, ERROR, IOFreq, false, BODY);
    
	SYS_all = CatStructFields(2, SYS_ref2_body, SYS_err2_body);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_all;
end

%% Group Data
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N{1,3}
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

%% Figure
close all ; clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 8])
movegui(fig, 'center')

pI = 1;
clear ax h
ax(1) = subplot(4,1,1); cla ; hold on ; ylim([0 1])
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,pI,:)), ...
        '.-', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.25)
    plot(GRAND.fly_stats(1).mean.IOFv.mean(:,1), GRAND.fly_stats(1).mean.IOGain.mean(:,pI), ...
        '.-k', 'MarkerSize', 17, 'LineWidth', 2)
    
ax(2) = subplot(4,1,2); cla ; hold on ; ylim([-300 150])
    yline(0, '--');
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,pI,:))), ...
        '.-', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.25)
    plot(GRAND.fly_stats(1).mean.IOFv.mean(:,1), rad2deg(GRAND.fly_stats(1).circ_mean.IOPhaseDiff.circ_mean(:,pI)), ...
        '.-k', 'MarkerSize', 17, 'LineWidth', 2)
    
ax(3) = subplot(4,1,3); cla ; hold on ; ylim([0 1.5])
    yline(1, '--');
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOFRF_error(:,pI,:)), ...
        '.-', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.25)
    plot(GRAND.fly_stats(1).mean.IOFv.mean(:,1), GRAND.fly_stats(1).mean.IOFRF_error.mean(:,pI), ...
        '.-k', 'MarkerSize', 17, 'LineWidth', 2)
    
ax(4) = subplot(4,1,4); cla ; hold on ; ylim([0 1])
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,pI,:)), ...
        '-', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.25)
    plot(GRAND.fly_stats(1).mean.Fv.mean(:,1), GRAND.fly_stats(1).mean.Cohr.mean(:,pI), ...
        '-k', 'MarkerSize', 17, 'LineWidth', 2)
    xlabel('Frequency (hz)')
    
set(ax, 'Color', 'none', 'LineWidth', 1, 'XScale', 'log', 'XLim', [0.2 25])
linkaxes(ax, 'x')


%% SAVE
disp('Saving...')
savedir = fullfile(rootdir, 'data'); mkdir(savedir)
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end