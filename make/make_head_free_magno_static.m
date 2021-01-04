function [] = make_head_free_magno_static(rootdir)
%% make_head_free_magno_static:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_static_wave';
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_static_wave_head_fixed';

%% Setup Directories %%
root.base = rootdir; clear rootdir
root.body = fullfile(root.base,'tracked_body');
root.reg = fullfile(root.base,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');
head_test = exist(root.head, 'file');

if head_test == 7
    filename = 'Static_HeadFree';
else
	filename = 'Static_HeadFixed';
end

% Select files
% [D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);
[D,I,N,U,T,~,~,basename] = GetFileData(root.body,'*.mat',false);


%% Get Data %%
close all
clc

Fs = 100;
Fc = 40;
func_length = 20;
startI = round(5000*0.5);
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
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
    if head_test == 7
        data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
    else % no head data
        data.head.head_data.angle = 0*data.body.bAngles;
    end
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,~]    = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, startI, false, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    lwing = rad2deg(data.benifly.LWing);
    rwing = rad2deg(data.benifly.RWing);
    lwing = hampel(data.benifly.Time, lwing);
    rwing = hampel(data.benifly.Time, rwing);
	lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
    % Get head & body anlges
    body = data.body.bAngles;
	head = data.head.head_data.angle;
    
    % Interpolate so all signals have the same times    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = -interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}   = nan;
    DATA.body{n}        = singal_attributes(Body, tintrp, [], n_detrend);
    DATA.head{n}        = singal_attributes(Head, tintrp, []);
    DATA.error{n}       = nan;
    DATA.dwba{n}    	= singal_attributes(dWBA, tintrp, []);
    DATA.lwing{n}    	= singal_attributes(LWing, tintrp, []);
    DATA.rwing{n}       = singal_attributes(RWing, tintrp, []);
    
    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Body','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
end

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end