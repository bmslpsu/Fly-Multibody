function [] = make_head_free_magno_SOS_amp(rootdir)
%% MakeData_SOS_v1_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v1';
filename = 'SOS_HeadFree_amp';

%% Setup Directories %%
root.daq = rootdir; clear rootdir
root.body = fullfile(root.daq,'tracked_body');
root.reg = fullfile(root.daq,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

% Load function files
FUNC = cell(N.amp, 1);
FUNC{1} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v1\function\All_position_function_SOS_Fs_50_T_20_freq_0.45_2.55_4.6_6.7_8.75_10.85_12.9_15_amp_1_1_1_1_1_1_1_1.mat');
FUNC{2} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v1\function\All_position_function_SOS_Fs_50_T_20_freq_0.35_2.15_4_5.8_7.6_9.4_11.25_13.05_amp_2_2_2_2_2_2_2_2.mat');
FUNC{3} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_amp_v1\function\All_position_function_SOS_Fs_50_T_20_freq_0.25_1.9_3.55_5.2_6.9_8.55_10.2_11.85_amp_3_3_3_3_3_3_3_3.mat');

%% Get Data %%
close all
clc

Fs = 100;
Fc = 20;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N.amp);
DATA = [I , splitvars(table(num2cell(zeros(N.file,5))))]; % store saccade objects
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error','dwba',};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.daq,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...
                            [basename{n} '.csv'])); % load head & wing angles from Benifly
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length,trigger, true, [], false, false);
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
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Body - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
%     [R,P] = corr(LWing, RWing);
%     x = LWing - mean(LWing);
%     y = RWing - mean(RWing);
%     [r,lags] = xcorr(x, y, 'normalized');
%     tlags = lags*(1/Fs);
%     subplot(2,1,1); cla ; hold on
%         plot(tintrp, x, 'b')
%         plot(tintrp, y, 'r')
%     subplot(2,1,2); cla ; hold on
%         plot(tlags, r, 'k', 'LineWidth', 1)
%         
%     pause
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}   = singal_attributes(Reference, tintrp);
    DATA.body{n}        = singal_attributes(Body, tintrp, [], n_detrend);
    DATA.head{n}        = singal_attributes(Head, tintrp, []);
    DATA.error{n}       = singal_attributes(Error, tintrp, [], n_detrend);
    DATA.dwba{n}    	= singal_attributes(dWBA, tintrp, []);
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, []);
    DATA.rwing{n}       = singal_attributes(RWing, tintrp, []);

%     subplot(1,1,1); cla ; hold on ; title(I.amp(n))
%         plot(DATA.reference{n}.Fv, DATA.reference{n}.mag.position, 'k')
%         pause
    
    % Debug plot.
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
    
    if 0
        figure (200)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, LWing, 'r', 'LineWidth', 1)
            plot(tintrp, RWing, 'b', 'LineWidth', 1)
            plot(tintrp, dWBA, 'k', 'LineWidth', 1)
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
     
    IOFreq = FUNC{I.amp(n)}.All.Freq;
    %IOFreq = flipud(FUNC{I.amp(n)}.All.Freq);
    SYS_ref2_head_body  = frf(tintrp, Reference, IOFreq, false, Body, Head);
    %pause
    %close all
    SYS_ref2_wing       = frf(tintrp, Reference, IOFreq, false, LWing, RWing, dWBA);
    SYS_head2_body_wing = frf(tintrp, Head, IOFreq, false, Body, dWBA);
    SYS_wing2_body      = frf(tintrp, dWBA, IOFreq, false, Body);
    SYS_left2_right  	= frf(tintrp, LWing, IOFreq, false, -RWing);
    
    SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_ref2_wing, ...
                                    SYS_head2_body_wing, SYS_wing2_body, SYS_left2_right);
    
    ALL{I.fly(n),I.amp(n)}(end+1,1) = SYS_all;
end

%% Group Data
clc
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N.amp
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
                
                GRAND.fly_stats(v).(stat_fields{s}).(fields{f}) = ...
                    system_stats(GRAND.fly_all(v).(stat_fields{s}).(fields{f}),3);
            end
        end
    end
    GRAND.all_trial(v) = structfun(@(x) system_stats(x,3), GRAND.all(v), 'UniformOutput', false);
end

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'ALL', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end