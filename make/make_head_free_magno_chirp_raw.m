function [] = make_head_free_magno_chirp_raw(rootdir)
%% make_head_free_magno_chirp:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_chirp_15';
exp_name = textscan(char(rootdir), '%s', 'delimiter', '_');
exp_typ = exp_name{1}{end-1}; % type of stimuli (vel or pos)
exp_ver = exp_name{1}{end}; % version of experiment (v1, v2, ...)

% clss = 'position';
clss = 'velocity';
filename = ['Chirp_HeadFree_' exp_typ '_' exp_ver '_' num2str(clss)];

%% Setup Directories %%
root.base = rootdir;
root.reg = fullfile(root.base,'registered');
root.body = fullfile(root.base,'tracked_body');
root.head = fullfile(root.reg,'tracked_head_tip');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);

%% Get Data
close all
clc

Fs = 100;
Fc = 40;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
[b,a] = butter(3, Fc/(Fs/2),'low');
DATA = [D , splitvars(table(num2cell(zeros(N.file,3))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','head'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.base,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
 	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles

    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT] 	= sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                    trigger, true, [], false, false);
    if length(TRIG.time_sync) == length(data.body.bAngles)
        % pass
    else
        [TRIG,PAT] = sync_pattern_trigger(daq_time, daq_pattern, func_length, ...
                        trigger, true, [], true, false); 
    end
 	trig_time = TRIG.time_sync;
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.bAngles;
    head = data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    %[b_pat, a_pat] = butter(3, 20 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = Reference - mean(Reference);
    Reference = 3.75*round(Reference/3.75);
    %Reference = filtfilt(b_pat, a_pat, Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    
    DATA.reference{n}       = Reference;
    DATA.body{n}            = Body;
    DATA.head{n}            = Head;
end

%% Time
pI = 1;
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 3])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,1);
Fv = GRAND.fly_stats.mean(v).Fv.mean;
cc = [0 0.1 0.7; 1 0 0; 1 0 0];
fly_lw = 0.25;
time = GRAND.fly_stats.mean(v).Time.mean;
ref = 15*chirp(time,0.1,20,6.5, 'logarithmic');
ref_vel = central_diff(ref, 1/Fs);
ax(1,1) = subplot(1,1,1); hold on
    plot(time, ref_vel, 'Color', [0 0 0], 'LineWidth', 1)
    %plot(time, squeeze(GRAND.fly_all.mean(1).State(:,pI,:)), ...
        %'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
%     plot(time, squeeze(GRAND.all(1).State(:,pI,:)), ...
%         'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
    [h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.fly_stats.mean(v).State.mean(:,pI),...
              GRAND.fly_stats.mean(v).State.std(:,pI), time, 1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
%     [h.patch(2,v),h.line(2,v)] = PlotPatch(GRAND.fly_stats.mean(v).State.mean(:,3),...
%               GRAND.fly_stats.mean(v).State.std(:,3), time, 1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
ax(1,1).Position(4) = 0.8*ax(1,1).Position(4);
ax(1,1).Position(2) = 1.5*ax(1,1).Position(2);
xlabel('Time (s)')
ylabel('Body (°)')
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0 20], 'YLim', 600*[-1 1], ...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')         
          
ax_freq = axes;
set(ax_freq, 'Color', 'none', 'YColor', 'none', 'XColor', 'k', 'Position', ax(1).Position, ...
    'XAxisLocation', 'top', 'XLim', [0.1 6.5], 'XScale', 'log', 'XTick', [0.1 1:6], 'FontSize', 8, ...
    'LineWidth', ax(1).LineWidth)
xlabel('Frequency (Hz)')

%% SAVE
time = tintrp;
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'time', '-v7.3')
disp('SAVING DONE')
end