function [] = make_head_free_magno_SOS_body_fixed(rootdir)
%% make_head_free_magno_SOS_body_fixed:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2_body_fixed';
exp_name = textscan(char(rootdir), '%s', 'delimiter', '_');
exp_typ = exp_name{1}{end-3}; % version of experiment (v1, v2, ...)
exp_ver = exp_name{1}{end-2}; % % type of stimuli (vel or pos)

clss = 'position';
% clss = 'velocity';
filename = ['SOS_BodyFixed_' exp_typ '_' exp_ver '_' num2str(clss)];

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.benifly = fullfile(root.base ,'tracked_head_wing');
root.head = fullfile(root.base ,'tracked_head_tip');
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
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

%% Get Data %%
close all
clc

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
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
 	%data.benifly = ImportBenifly(fullfile(root.benifly, ...  % load head & wing angles from Benifly
                            %[basename{n} '.csv']));
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, [], false, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    %lwing = rad2deg(data.benifly.LWing);
    %rwing = rad2deg(data.benifly.RWing);
    %lwing = hampel(data.benifly.Time, lwing);
    %%rwing = hampel(data.benifly.Time, rwing);
	%lwing = filtfilt(b,a,lwing);
    %rwing = filtfilt(b,a,rwing);
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
	head = -data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    %Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    Reference = 3.75*(round(Reference/3.75));
    
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Head;
    %LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    %RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    %dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    DATA.body_saccade{n}    = [];
    DATA.reference{n}       = singal_attributes(Reference, tintrp);
    DATA.body{n}            = [];
    DATA.head{n}            = singal_attributes(Head, tintrp);
    DATA.error{n}           = singal_attributes(Error, tintrp);
    %DATA.dwba{n}            = singal_attributes(dWBA, tintrp, [], n_detrend);
    %DATA.dwba{n}            = singal_attributes(DATA.dwba{n}.detrend, tintrp, [], []);
    %DATA.lwing{n}           = singal_attributes(LWing, tintrp, []);
    %DATA.rwing{n}           = singal_attributes(RWing, tintrp, []);
        
    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            %plot(tintrp, DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Reference','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = sort(FUNC{I{n,3}}.All.Freq, 'ascend');
    REF = DATA.reference{n}.(clss);
    HEAD = DATA.head{n}.(clss);
    %dWBA = DATA.dwba{n}.(clss);
    %LWING = DATA.lwing{n}.(clss);
    %RWING = DATA.rwing{n}.(clss);
    
    SYS_ref2_head = frf(tintrp, REF , IOFreq, false, HEAD);
    %SYS_ref2_wing = frf(tintrp, REF, IOFreq, false, dWBA);
    %SYS_wing2_body = frf(tintrp, dWBA, IOFreq, false, BODY);
    
	%SYS_all = CatStructFields(2, SYS_ref2_head);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_ref2_head;
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

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end