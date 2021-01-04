function [] = make_head_free_magno_SOS_all_head_fixed(rootdir)
%% make_head_free_magno_SOS_all_head_fixed:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2_head_fixed';
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v3_head_fixed';
exp_name = textscan(char(rootdir), '%s', 'delimiter', '_');
exp_typ = exp_name{1}{end-3}; % version of experiment (v1, v2, ...)
exp_ver = exp_name{1}{end-2}; % % type of stimuli (vel or pos)

% clss = 'position';
clss = 'velocity';
filename = ['SOS_HeadFixed_' exp_typ '_' exp_ver '_' num2str(clss)];

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.body = fullfile(root.base,'tracked_body');
root.reg = fullfile(root.base,'registered');
root.benifly = fullfile(root.reg ,'tracked_wing');
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

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.body,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

%% Get Data %%
close all
clc

Fs = 100;
Fc = 20;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N{1,3});
DATA = [I , splitvars(table(num2cell(zeros(N.file,7))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error','dwba','lwing','rwing'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
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
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Error   = Reference - Body;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}   = singal_attributes(Reference, tintrp);
    DATA.body{n}        = singal_attributes(Body, tintrp, [], n_detrend);
    DATA.error{n}       = singal_attributes(Error, tintrp, [], n_detrend);
    DATA.dwba{n}    	= singal_attributes(dWBA, tintrp, []);
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, []);
    DATA.rwing{n}       = singal_attributes(RWing, tintrp, []);
    
    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            plot(tintrp, DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Reference','Body','\DeltaWBA', 'Orientation', 'horizontal');
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
    dWBA = DATA.dwba{n}.(clss);
    LWING = DATA.lwing{n}.(clss);
    RWING = DATA.rwing{n}.(clss);
    
    SYS_ref2_body = frf(tintrp, REF , IOFreq, false, BODY);
    SYS_ref2_wing = frf(tintrp, REF, IOFreq, false, LWING, RWING, dWBA);
    SYS_wing2_body = frf(tintrp, dWBA, IOFreq, false, BODY);
    %SYS_left2_right = frf(tintrp, LWING, IOFreq, false, RWING);
    
    SYS_all = CatStructFields(2, SYS_ref2_body, SYS_ref2_wing, SYS_wing2_body);
    
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

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end