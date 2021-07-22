function [] = make_head_free_magno_chirp(rootdir)
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
filename = ['SS_HeadFree_' exp_typ '_' exp_ver '_' num2str(clss)];

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

% Body saccade detection parameters
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
ALL = cell(N.fly,N{1,4});
DATA = [D , splitvars(table(num2cell(zeros(N.file,8))))];
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error',...
    'dwba','lwing','rwing','body_saccade'};
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
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Body - Head;
    
    % Detect & remove saccades
    body_scd = saccade_v1(Body, tintrp, scd.thresh, scd.true_thresh, scd.Fc_detect, ...
                            scd.Fc_ss, scd.amp_cut, scd.dur_cut , scd.direction, scd.pks, ...
                            scd.sacd_length, scd.min_pkdist, scd.min_pkwidth, scd.min_pkprom, ...
                            scd.min_pkthresh, scd.boundThresh, false);
%     figure (1)
%     pause
%     close all
    
    % Store signals
    n_detrend = 4;
    DATA.body_saccade{n}    = body_scd;
    DATA.reference{n}       = singal_attributes(Reference, tintrp);
    DATA.body{n}            = singal_attributes(body_scd.shift.IntrpPosition, tintrp, [], n_detrend);
    DATA.head{n}            = singal_attributes(Head, tintrp);
    DATA.error{n}           = singal_attributes(Error, tintrp);

    % Debug plot
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            leg = legend('Reference','Body', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    REF = DATA.reference{n}.(clss);
    BODY = DATA.body{n}.(clss);
    %HEAD = DATA.head{n}.(clss);
    
    bin_sz = 5;
    SYS_ref2_body_head = frf_chirp(tintrp, REF, bin_sz, false, BODY);

    %SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_head2_body, SYS_ref2_wing, SYS_wing2_body);
    
    ALL{I.fly(n),1}(end+1,1) = SYS_ref2_body_head;
end

%% Group Data
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N.A
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

%%
Fv = GRAND.fly_stats.mean(v).Fv.mean;
plot(Fv, squeeze(GRAND.all(1).Mag(:,pI,:)), ...
        'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
xlim([0 6.5])

%% BODE
n_cond = N.A;
pI = 1;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond 5*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_cond);
Fv = GRAND.fly_stats.mean(v).Fv.mean;
cc = [0 0.1 0.7; 1 0 0; 1 0 0];
fly_lw = 0.25;
for v = 1
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(5,n_cond,subI(1)); hold on
        plot(Fv, squeeze(GRAND.fly_all.mean(1).Mag(:,pI,:)), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(GRAND.fly_stats.mean(v).Mag.mean(:,pI),...
                  GRAND.fly_stats.mean(v).Mag.std(:,pI), Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(2,v) = subplot(5,n_cond,subI(2)); hold on
        plot(Fv, squeeze(GRAND.fly_all.mean(1).Gain(:,pI,:)), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(GRAND.fly_stats.mean(v).Gain.mean(:,pI),...
                  GRAND.fly_stats.mean(v).Gain.std(:,pI), Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(3,v) = subplot(5,n_cond,subI(3)); hold on
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        plot(Fv, rad2deg(squeeze(GRAND.fly_all.mean(1).PhaseDiff(:,pI,:))), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(rad2deg(GRAND.fly_stats.mean(v).PhaseDiff.mean(:,pI)),...
                  rad2deg(GRAND.fly_stats.mean(v).PhaseDiff.std(:,pI)), Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(4,v) = subplot(5,n_cond,subI(4)); hold on
        yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        plot(Fv, squeeze(GRAND.fly_all.mean(1).FRF_error(:,pI,:)), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        [h.patch(4,v,1),h.line(4,v,1)] = PlotPatch(GRAND.fly_stats.mean(v).FRF_error.mean(:,pI),...
                  GRAND.fly_stats.mean(v).FRF_error.std(:,pI), Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(5,v) = subplot(5,n_cond,subI(5)); hold on
        plot(Fv, squeeze(GRAND.fly_all.mean(1).Cohr(:,pI,:)), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        [h.patch(5,v,1),h.line(5,v,1)] = PlotPatch(GRAND.fly_stats.mean(v).Cohr.mean(:,pI),...
                  GRAND.fly_stats.mean(v).Cohr.std(:,pI), Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);

end
linkaxes(ax, 'x')

% delete(h.patch)

set(h.line, 'Marker', 'none','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.1 6.5],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1:6])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Compensation error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 6])
set(ax(2,1:end-1),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-150 150])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')
% 
% set(ax,'XScale','log')

align_Ylabels(fig)

%% BODE Grand
tt = GRAND.fly_stats.mean(v).Time.mean;
REF = GRAND.fly_stats.mean(v).refState.mean;
BODY = GRAND.fly_stats.mean(v).State.mean(:,1);
HEAD = GRAND.fly_stats.mean(v).State.mean(:,2);
test = frf_chirp(tt, REF, bin_sz, false, BODY, HEAD);

n_cond = N.A;
pI = 1;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond 5*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_cond);
Fv = GRAND.fly_stats.mean(v).Fv.mean;
% cc = [0 0.1 0.7];
fly_lw = 0.25;
for v = 1
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(5,n_cond,subI(1)); hold on
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(GRAND.fly_stats.mean(v).Mag.mean(:,pI),...
                  GRAND.fly_stats.mean(v).Mag.std(:,pI)*0, Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(2,v) = subplot(5,n_cond,subI(2)); hold on
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(squeeze(test.Gain(:,pI,:)),...
                  GRAND.fly_stats.mean(v).Gain.std(:,pI)*0, Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(3,v) = subplot(5,n_cond,subI(3)); hold on
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(rad2deg(squeeze(test.PhaseDiff(:,pI,:))),...
                  rad2deg(squeeze(test.Gain(:,pI,:)))*0, Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(4,v) = subplot(5,n_cond,subI(4)); hold on
        yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(4,v,1),h.line(4,v,1)] = PlotPatch(squeeze(test.FRF_error(:,pI,:)),...
                  squeeze(test.FRF_error(:,pI,:))*0, Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
              
    ax(5,v) = subplot(5,n_cond,subI(5)); hold on
        [h.patch(5,v,1),h.line(5,v,1)] = PlotPatch(squeeze(test.Cohr(:,pI,:)),...
                  squeeze(test.Cohr(:,pI,:))*0, Fv, 1, 1, cc(pI,:), 0.7*cc(pI,:), 0.2, 1);
end
linkaxes(ax, 'x')

% delete(h.patch)

set(h.line, 'Marker', 'none','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.1 6.5],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1:6])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Compensation error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 6])
set(ax(2,1:end-1),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-150 150])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')
% 
% set(ax,'XScale','log')

align_Ylabels(fig)

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end