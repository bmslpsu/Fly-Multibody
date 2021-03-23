function [] = Saccade_window()
%% Saccade_window:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'DATA','D','I','U','N')

%% Get saccades
clearvars -except DATA N U

SACCADE.body_pos = cell(N.freq,1);
SACCADE.body_vel = cell(N.freq,1);
SACCADE.head_pos = cell(N.freq,1);
SACCADE.head_vel = cell(N.freq,1);
SACCADE.dwba_pos = cell(N.freq,1);
SACCADE.dwba_vel = cell(N.freq,1);
freqI = findgroups(DATA.freq);
scd_win = 0.2;
norm_dir = false;
for n = 1:N.file
    scd = DATA.body_saccade{n};
    if scd.count > 0
        % Body
        [scds,~,scd_time,~] = getSaccade(scd, scd.position, scd_win, 0.5, norm_dir);
        scd_time = scd_time{1};
        scds = cellfun(@(x) x - median(x(1:20)), scds, 'UniformOutput', false);
        scds = cat(2,scds{:});
        sI = size(SACCADE.body_pos{freqI(n)},2) + 1;
        winI = sI + size(scds,2) - 1;
        SACCADE.body_pos{freqI(n)}(:,sI:winI) = scds;
        
        [scds,~,~,~] = getSaccade(scd, scd.velocity, scd_win, scd_win, norm_dir);
        scds = cat(2,scds{:});
        SACCADE.body_vel{freqI(n)}(:,sI:winI) = scds;
        
        % Head
        [scds,~,~,~] = getSaccade(scd, DATA.head{n}.position, scd_win, scd_win, norm_dir);
        scds = cellfun(@(x) x - median(x), scds, 'UniformOutput', false);
        scds = cat(2,scds{:});
        SACCADE.head_pos{freqI(n)}(:,sI:winI) = scds;
        
        [scds,~,~,~] = getSaccade(scd, DATA.head{n}.velocity, scd_win, scd_win, norm_dir);
        scds = cat(2,scds{:});
        SACCADE.head_vel{freqI(n)}(:,sI:winI) = scds;
        
        % dWBA
        [scds,~,~,~] = getSaccade(scd, DATA.dwba{n}.position, scd_win, scd_win, norm_dir);
        scds = cellfun(@(x) x - median(x), scds, 'UniformOutput', false);
        scds = cat(2,scds{:});
        SACCADE.dwba_pos{freqI(n)}(:,sI:winI) = scds;
        
        [scds,~,~,~] = getSaccade(scd, DATA.dwba{n}.velocity, scd_win, scd_win, norm_dir);
        scds = cat(2,scds{:});
        SACCADE.dwba_vel{freqI(n)}(:,sI:winI) = scds;
    end 
end
SACCADE_all = structfun(@(x) cat(2,x{:}), SACCADE, 'UniformOutput', false);
SACCADE.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    SACCADE, 'UniformOutput', false);
SACCADE_all.stats = structfun(@(x) basic_stats(x,2), SACCADE_all, 'UniformOutput', false);

%% Plot saccades all
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 4])
movegui(fig, 'center')
clear ax H
ax = gobjects(2,1);
grand_lw = 1;
body_color = [1 0 0];
head_color = [0 0 1];
wing_color = [0 1 0];
n_std = 1;
ax(1) = subplot(2,1,1); hold on ; cla
%     [~] = PlotPatch(SACCADE_all.stats.body_pos.mean, SACCADE_all.stats.body_pos.std, ...
%         scd_time, n_std, 1, body_color, body_color, 0.3, grand_lw);
    [~] = PlotPatch(SACCADE_all.stats.head_pos.mean, SACCADE_all.stats.head_pos.std, ...
        scd_time, n_std, 1, head_color, head_color, 0.3, grand_lw);
    [~] = PlotPatch(SACCADE_all.stats.dwba_pos.mean, SACCADE_all.stats.dwba_pos.std, ...
        scd_time, n_std, 1, wing_color, wing_color, 0.3, grand_lw);

ax(2) = subplot(2,1,2); hold on ; cla
%     [~] = PlotPatch(SACCADE_all.stats.body_vel.mean, SACCADE_all.stats.body_vel.std, ...
%         scd_time, n_std, 1, body_color, body_color, 0.3, grand_lw);
    [~] = PlotPatch(SACCADE_all.stats.head_vel.mean, SACCADE_all.stats.head_vel.std, ...
        scd_time, n_std, 1, head_color, head_color, 0.3, grand_lw);
    [~] = PlotPatch(SACCADE_all.stats.dwba_vel.mean, SACCADE_all.stats.dwba_vel.std, ...
        scd_time, n_std, 1, wing_color, wing_color, 0.3, grand_lw);
    
set(ax, 'XLim', 0.2*[-1 1])
% set(ax(1,:), 'YLim', [-10 80])
% set(ax(2,:), 'YLim', [-200 1200])
% set(ax(3,:), 'YLim', 40000*[-1 1])
set(ax, 'LineWidth', 1, 'Color', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
% YLabelHC = get(ax(3,1), 'YLabel');
% set([YLabelHC], 'String', 'Accleratin (°/s^{2})')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Plot saccades
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 8 4.5])
movegui(fig, 'center')
clear ax H
ax = gobjects(2,N.freq);
fly_lw = 0.5;
grand_lw = 1;
body_color = [1 0 0];
head_color = [0 0 1];
wing_color = [0 1 0];
n_std = 1;
for f = 1:N.freq
    rowI = f + (0:2)*N.freq;
    ax(1,f) = subplot(2,N.freq,rowI(1)); hold on ; cla
        [~] = PlotPatch(SACCADE.fly_stats.body_pos(f).mean, SACCADE.fly_stats.body_pos(f).std, ...
            scd_time, n_std, 1, body_color, body_color, 0.3, grand_lw);
        [~] = PlotPatch(SACCADE.fly_stats.head_pos(f).mean, SACCADE.fly_stats.head_pos(f).std, ...
            scd_time, n_std, 1, head_color, head_color, 0.3, grand_lw);
  
    ax(2,f) = subplot(2,N.freq,rowI(2)); hold on ; cla
        [~] = PlotPatch(SACCADE.fly_stats.body_vel(f).mean, SACCADE.fly_stats.body_vel(f).std, ...
            scd_time, n_std, 1, body_color, body_color, 0.3, grand_lw);
        [~] = PlotPatch(SACCADE.fly_stats.head_vel(f).mean, SACCADE.fly_stats.head_vel(f).std, ...
            scd_time, n_std, 1, head_color, head_color, 0.3, grand_lw);
%         plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_vel{v}, ...
%             'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw) 
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
% linkaxes(ax(3,:), 'y')
set(ax, 'XLim', 0.2*[-1 1])
set(ax(1,:), 'YLim', 80*[-1 1])
set(ax(2,:), 'YLim', [-200 1500])
% set(ax(3,:), 'YLim', 40000*[-1 1])
set(ax(1:2,:), 'XTickLabel', [])
set(ax(1:2,:), 'XColor', 'none')
set(ax(:,2:N.freq), 'YTickLabel', [])
set(ax(:,2:N.freq), 'YColor', 'none')
% set(ax(3,2:N.freq), 'XColor', 'none')
set(ax, 'LineWidth', 1, 'Color', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
% YLabelHC = get(ax(3,1), 'YLabel');
% set([YLabelHC], 'String', 'Accleratin (°/s^{2})')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Save
% savedir = 'E:\DATA\Magno_Data\Multibody\processed';
% filename = ['Saccade_win_' clss '.mat'];
% save(fullfile(savedir, filename), 'All','U','N','-v7.3')

end