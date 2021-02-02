function [] = make_head_free_magno_SOS_body_fixed(rootdir)
%% make_head_free_magno_SOS_body_fixed: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_body_fixed';
filename = 'SOS_HeadFree_body_fixed';

%% Setup Directories %%
root.daq = rootdir; clear rootdir
root.benifly = fullfile(root.daq ,'tracked_head_wing');
root.head = fullfile(root.daq ,'tracked_head_tip');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

% Load function files
FUNC = cell(N{1,3}, 1);
FUNC{1} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_body_fixed\functions\ALL_position_function_SOS_Fs_50_T_20_freq_1_3.1_5.3_7.4_9.6_amp_9_5_3_2_1.mat');

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
DATA = [I , splitvars(table(num2cell(zeros(N.file,5))))]; % store saccade objects
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error','dwba',};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.daq,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
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
    lwing = rad2deg(-data.benifly.LWing);
    rwing = rad2deg(-data.benifly.RWing);
    lwing = hampel(data.benifly.Time, lwing);
    rwing = hampel(data.benifly.Time, rwing);
	lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
    % Get pattern, head, & body angles
    pat = 3.75*PAT.pos;
	head = -data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference  - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}   = singal_attributes(Reference, tintrp);
    DATA.head{n}        = singal_attributes(Head, tintrp, []);
    DATA.error{n}       = singal_attributes(Error, tintrp, [], n_detrend);
    DATA.dwba{n}    	= singal_attributes(dWBA, tintrp, []);
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, []);
    DATA.rwing{n}       = singal_attributes(RWing, tintrp, []);

    % Debug plot.
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
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
     
    IOFreq = FUNC{I{n,3}}.All.Freq;
    %IOFreq = flipud(FUNC{I.amp(n)}.All.Freq);
    SYS_ref2_head  = frf(tintrp, Reference, IOFreq, false, Head);
    %pause
    %close all
    SYS_ref2_wing       = frf(tintrp, Reference, IOFreq, false, LWing, RWing, dWBA);
    SYS_head2_wing      = frf(tintrp, Head, IOFreq, false, dWBA);
    SYS_left2_right  	= frf(tintrp, LWing, IOFreq, false, -RWing);
    
    SYS_all = CatStructFields(2, SYS_ref2_head, SYS_ref2_wing, ...
                                    SYS_head2_wing, SYS_left2_right);
    
    ALL{I.fly(n),I{n,3}}(end+1,1) = SYS_all;
end

%% Group Data
clc
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