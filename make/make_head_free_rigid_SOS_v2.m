function [] = make_head_free_rigid_SOS_v2(rootdir)
%% make_head_free_magno_SOS_amp:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_SOS_v2';
ver = textscan(char(rootdir), '%s', 'delimiter', '_');
ver = ver{1}{end}; % version of experiment (v1, v2, ...)

clss = 'position';
filename = ['SOS_HeadFree_amp_' ver '_' num2str(clss) '_rigid'];

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.filtvid = fullfile(root.base ,'filtvid');
root.benifly = fullfile(root.filtvid ,'tracked_head_wing');
root.head = fullfile(root.base ,'tracked_head');

% Select files
% [D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
[D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

%% Get Data %%
close all
clc

Fs = 100;
Fc = 20;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,1);
DATA = [I , splitvars(table(num2cell(zeros(N.file,5))))]; % store saccade objects
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error','dwba',};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'hAngles'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...  % load head & wing angles from Benifly
                            [basename{n} '.csv']));
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, [], true, false);
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
	head = data.head.hAngles;
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}   = singal_attributes(Reference, tintrp);
    DATA.head{n}        = singal_attributes(Head, tintrp, []);
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
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, 5*DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Reference','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
    REF = DATA.reference{n}.(clss);
    HEAD = DATA.head{n}.(clss);
    dWBA = DATA.dwba{n}.(clss);
    LWING = DATA.lwing{n}.(clss);
    RWING = DATA.rwing{n}.(clss);
    
    SYS_ref2_head = frf(tintrp, REF , IOFreq, false, HEAD);
    SYS_ref2_wing = frf(tintrp, REF, IOFreq, false, dWBA);
    SYS_head2_wing = frf(tintrp, HEAD, IOFreq, false, dWBA);
    
    SYS_all = CatStructFields(2, SYS_ref2_head, SYS_ref2_wing, SYS_head2_wing);
    
    ALL{I.fly(n),1}(end+1,1) = SYS_all;
end

%% Group Data
% clc
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1
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
savedir = 'H:\DATA';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
     'DATA', 'ALL', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end