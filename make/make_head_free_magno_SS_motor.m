function [] = make_head_free_magno_SS_motor(rootdir)
%% make_head_free_magno_SS_motor:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

warning('off', 'signal:findpeaks:largeMinPeakHeight')

clss = 'position';
% clss = 'velocity';

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250_motor_drum';

[~,exp_name,~] = fileparts(rootdir);
exp_name = textscan(char(exp_name), '%s', 'delimiter', '_');
exp_name = string(exp_name{1});
exp_name = char(strjoin(exp_name(2:end), '_'));
filename = [exp_name '_' num2str(clss)];

%% Setup Directories
root.base = rootdir;
% root.body = fullfile(root.base,'tracked_body');
root.reg = fullfile(root.base,'registered');
root.body = fullfile(root.reg,'reg_angles');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');
root.func = fullfile(root.base ,'function');
root.replay = fullfile(root.base ,'replay');

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

% Load replay file
func_list = dir(root.replay);
func_list = func_list(~[func_list.isdir]);
n_cond = length(func_list);
Replay = cell(n_cond,1);
for f = 1:n_cond
    Replay{f} = load(fullfile(root.replay, func_list(f).name));
    Replay{f}.name = func_list(f).name;
end
Replay = cellfun(@(x) x.replay, Replay, 'UniformOutput', false);

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);

%% Get Data
close all
clc

Fs = 100;
func_length = 10;
tintrp = (0:(1/Fs):func_length)';
debug = false;
ALL = cell(N.fly,N{1,3});
DATA = [D , splitvars(table(num2cell(zeros(N.file,8))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error',...
    'dwba','lwing','rwing','body_saccade'};
Extra = [];
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p','t_v'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'angles'); % load body angles
    % data.head = load(fullfile(root.head, [basename{n} '.mat']),'hAngles'); % load head angles
    data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
    % data.benifly = ImportBenifly(fullfile(root.benifly, ...  % load head & wing angles from Benifly
    %                         [basename{n} '.csv']));
    
    % Get head data & filter
    Fc_low = D.freq(n)*1.5;
    [b_low, a_low] = butter(5, Fc_low / (Fs/2), 'low');
    Fc_high = 0.3;
    [b_high, a_high] = butter(5, Fc_high / (Fs/2), 'high');
    head = hampel(data.daq.t_v, data.head.head_data.angle);
    head = filtfilt(b_low, a_low, head);
    head = filtfilt(b_high, a_high, head);
    
    % Get body data
    body = data.body.angles;
    
    % Put in DAQ reference frame
    thresh = FUNC{I.freq(n)}.All.Amp*0.5;
    %[first_peak] = find_first_peak(body, data.daq.t_v, 15, thresh, false);
    first_peak = 55;
    
    daq_time = data.daq.t_p;
    trigger = round(data.daq.data(:,1));
    [~,locs] = findpeaks(trigger, 'MinpeakHeight', 1);
    trig_time_raw = daq_time(locs);
    [~,first_peak_daq] = min(abs(trig_time_raw(first_peak) - daq_time));
    
    % Get synced frame times and pattern data
    daq_pattern = data.daq.data(:,2);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, first_peak_daq - 800, false, false);

    trig_time = TRIG.time_sync(1:end);
    
    % % Filter wing angles
    % lwing = rad2deg(data.benifly.LWing);
    % rwing = rad2deg(data.benifly.RWing);
    % lwing = hampel(data.benifly.Time, lwing);
    % rwing = hampel(data.benifly.Time, rwing);
    % lwing = filtfilt(b,a,lwing);
    % rwing = filtfilt(b,a,rwing);
    
    % Get pattern, head, & body anlges
    % pat = 3.75*PAT.pos;
    
    % % Interpolate so all signals have the same times
    % [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    % Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    % %Reference = filtfilt(b_pat, a_pat, Reference);
    % Reference = Reference - mean(Reference);
    % Reference = 3.75*(round(Reference/3.75));
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    % LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    % RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    % dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    % DATA.reference{n}       = singal_attributes(Reference, tintrp);
    DATA.body{n}            = singal_attributes(Body, tintrp, 20, 0);
    DATA.head{n}            = singal_attributes(Head, tintrp, 20, 0);
    % DATA.dwba{n}            = singal_attributes(dWBA, tintrp, 20, 0);
    % DATA.lwing{n}           = singal_attributes(LWing, tintrp, 20);
    % DATA.rwing{n}           = singal_attributes(RWing, tintrp, 20);
    
    %Error                   = -DATA.body{n}.position - DATA.head{n}.position;
	%DATA.error{n}           = singal_attributes(Error, tintrp);

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
    
    IOFreq = D.freq(n);
    %REF = DATA.reference{n}.(clss);
    BODY = DATA.body{n}.(clss);
    HEAD = DATA.head{n}.(clss);
    %ERROR = DATA.error{n}.(clss);
    %dWBA = DATA.dwba{n}.position; % rerun for this
    %LWING = DATA.lwing{n}.(clss);
    %RWING = DATA.rwing{n}.(clss);

    %SYS_ref2_head_body = frf(tintrp, REF , IOFreq, false, BODY, HEAD);
	SYS_body2_head = frf(tintrp, BODY, IOFreq, false, HEAD);
    %SYS_ref2_wing = frf(tintrp, DATA.reference{n}.position, IOFreq, false, dWBA);
    %SYS_wing2_body = frf(tintrp, dWBA, IOFreq, false, BODY);
 	%SYS_err2_head_body = frf(tintrp, ERROR, IOFreq, false, BODY, HEAD);
    
	%SYS_all = CatStructFields(2, SYS_body2_head, SYS_wing2_body);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_body2_head;
    Extra.Body_Freq(:,n) = DATA.body{n}.mag.velocity;
    Extra.Head_Freq(:,n) = DATA.head{n}.mag.velocity;
    
%     subplot(3,1,1) ; cla ; hold on ; title(D.freq(n))
%         plot(tintrp, Body, 'r')
%         plot(tintrp, Head, 'b')
%     subplot(3,1,2) ; cla ; hold on
%         plot(DATA.body{n}.Fv, DATA.body{n}.mag.position, 'r')
%         plot(DATA.body{n}.Fv, DATA.head{n}.mag.position, 'b')
%         xlim([0 20])
%     subplot(3,1,3) ; cla ; hold on
%         plot(DATA.body{n}.Fv, SYS_body2_head.Cohr, 'k')
%         xlim([0 20])
%         ylim([0 1])
%         
%  	pause
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

%% Check
m = 3;
clf
subplot(4,1,1) ; cla ; hold on
    body_all = squeeze(GRAND.all(m).refState(:,1,:));
    plot(tintrp, body_all - mean(body_all(1,:)), 'r', 'Color', [1 0 0 0.4], 'LineWidth', 0.75)
    plot(Replay{1}.time, Replay{1}.pos.body_sine(:,m) - Replay{1}.pos.body_sine(1,m), 'k', 'LineWidth', 0.5)
    
subplot(4,1,2) ; cla ; hold on
    plot(DATA.body{n}.Fv, Extra.Body_Freq, 'r', 'Color', [1 0 0 0.4], 'LineWidth', 0.75)
    plot(Replay{1}.Fv, Replay{1}.freq.vel.body_sine.mag(:,m), 'k', 'LineWidth', 0.5)
    xlim([0.25 15])
    set(gca, 'XScale', 'log')

subplot(4,1,3) ; cla ; hold on
    plot(squeeze(GRAND.all(m).Time(:,1,:)), squeeze(GRAND.all(m).State(:,1,:)), ...
        'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
    plot(GRAND.all_trial(m).Time.mean(:,1), GRAND.all_trial(m).State.mean(:,1), ...
        'Color', [0 0 1 1], 'LineWidth', 0.75)
    ylim(10*[-1 1])
    
subplot(4,1,4) ; cla ; hold on

%     plot(DATA.head{n}.Fv, Extra.Head_Freq, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
%     plot(GRAND.all_trial(1).Fv.mean(:,1), GRAND.all_trial(1).Mag.mean(:,1), ...
%         'Color', [0 0 1 1], 'LineWidth', 0.75)


    plot(squeeze(GRAND.all(m).Fv(:,1,:)), squeeze(GRAND.all(m).Mag(:,1,:)), ...
        'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
    plot(GRAND.all_trial(m).Fv.mean(:,1), GRAND.all_trial(m).Mag.mean(:,1), ...
        'Color', [0 0 1 1], 'LineWidth', 0.75)
    xlim([0.25 15])
    set(gca, 'XScale', 'log')

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'Replay', 'DATA', 'GRAND', 'FLY', 'Replay', 'Extra', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end