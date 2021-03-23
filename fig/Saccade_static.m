function [] = Saccade_static()
%% Saccade_static: compare saccade rate and dynamics of data sets
warning('off', 'signal:findpeaks:largeMinPeakHeight')
root = 'E:\DATA\Magno_Data\Multibody';
[FILES,PATH] = uigetfile({'*.mat'}, 'Select data files', root, 'MultiSelect','on');
FILES = cellstr(FILES);
n_file = length(FILES);
ALL = cell(n_file,1);
labels = string(zeros(n_file,1));
for n = 1:n_file
    ALL{n} = load(fullfile(PATH,FILES{n}),'DATA','D','I','U','N');
    filedata = textscan(char(FILES{n}), '%s', 'delimiter', '_');
    filedata = filedata{1};
    ALL{n}.name = [filedata{1} '_' filedata{2}];
    labels(n) = ALL{n}.name;
end

%% Get saccades
clearvars -except ALL n_file FILES labels
clc

dataset = ALL{1};
val = [22.5 30 60];
valI = any(dataset.DATA.wave == val,2);
dataset.DATA = dataset.DATA(valI,:);
dataset.D = dataset.D(valI,:);

use_val = false;
time_win = 0.2;
name = 'body';
get_name = ["head", "dwba"];
[stats] = get_saccade_stats(dataset, name, use_val, time_win, get_name);
%%
use_val = false;
time_win = 0.15;
name = 'head';
get_name = 'body';
[stats] = get_saccade_stats(dataset, name, use_val, time_win, get_name);

%% Plot
if use_val
    n_val = size(stats.val_stats.head_position,2);
else
	n_val = 1;
end
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 2*[2 2 n_val*2 3*2])
movegui(fig, 'center')
clear ax H
ax = gobjects(3,n_val);
fly_lw = 0.25;
vel_lw = 1;
body_color = [0.9 0 0];
head_color = [0 0.4 1];
wing_color = [1 0 1];
n_std = 1;
fly_stop = nan;
for v = 1:n_val
    rowI = v + (0:2)*n_val;
    ax(1,v) = subplot(3,n_val,rowI(1)); hold on ; cla
        plot(stats.time, fly_stop*stats.fly_mean.body_position{v}, ...
            'Color', [0.7*body_color 0.5], 'LineWidth', fly_lw)
        
        plot(stats.time, fly_stop*stats.fly_mean.head_position{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)

        plot(stats.time, 1*fly_stop*stats.fly_mean.dwba_position{v}, ...
            'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(stats.val_stats.body_position(v).mean, stats.val_stats.body_position(v).std, ...
            stats.time, n_std, 1, body_color, body_color, 0.3, vel_lw);
        
        [~] = PlotPatch(stats.val_stats.head_position(v).mean, stats.val_stats.head_position(v).std, ...
            stats.time, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
        [~] = PlotPatch(1*stats.val_stats.dwba_position(v).mean, stats.val_stats.dwba_position(v).std, ...
            stats.time, n_std, 1, wing_color, wing_color, 0.3, vel_lw);
        
    ax(2,v) = subplot(3,n_val,rowI(2)); hold on ; cla
        plot(stats.time, fly_stop*stats.fly_mean.body_velocity{v}, ...
            'Color', [0.7*body_color 0.5], 'LineWidth', fly_lw)
        
        plot(stats.time, fly_stop*stats.fly_mean.head_velocity{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
        plot(stats.time, fly_stop*stats.fly_mean.dwba_velocity{v}, ...
            'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(stats.val_stats.body_velocity(v).mean, stats.val_stats.body_velocity(v).std, ...
            stats.time, n_std, 1, body_color, body_color, 0.3, vel_lw);
        
        [~] = PlotPatch(stats.val_stats.head_velocity(v).mean, stats.val_stats.head_velocity(v).std, ...
            stats.time, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
        [~] = PlotPatch(stats.val_stats.dwba_velocity(v).mean, stats.val_stats.dwba_velocity(v).std, ...
            stats.time, n_std, 1, wing_color, wing_color, 0.3, vel_lw);
        
    ax(3,v) = subplot(3,n_val,rowI(3)); hold on ; cla
        plot(stats.time, fly_stop*stats.fly_mean.body_acceleration{v}, ...
            'Color', [0.7*body_color 0.5], 'LineWidth', fly_lw)
        
        plot(stats.time, fly_stop*stats.fly_mean.dwba_acceleration{v}, ...
            'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)
        
        plot(stats.time, fly_stop*stats.fly_mean.head_acceleration{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(stats.val_stats.body_acceleration(v).mean, stats.val_stats.body_acceleration(v).std, ...
            stats.time, n_std, 1, body_color, body_color, 0.3, vel_lw);
        
        [~] = PlotPatch(stats.val_stats.head_acceleration(v).mean, stats.val_stats.head_acceleration(v).std, ...
            stats.time, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
        [~] = PlotPatch(stats.val_stats.dwba_acceleration(v).mean, stats.val_stats.dwba_acceleration(v).std, ...
            stats.time, n_std, 1, wing_color, wing_color, 0.3, vel_lw); 
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
set(ax, 'XLim', 0.15*[-1 1])
set(ax(1,:), 'YLim', [-8 50])
set(ax(2,:), 'YLim', [-300 1000])
set(ax(3,:), 'YLim', 30000*[-1 1])
set(ax(1:2,:), 'XTickLabel', [])
set(ax(1:2,:), 'XColor', 'none')
set(ax(:,2:n_val), 'YTickLabel', [])
set(ax(:,2:n_val), 'YColor', 'none')
set(ax(3,2:n_val), 'XColor', 'none')
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Accleratin (°/s^{2})')
XLabelHC = get(ax(3,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')


%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
dataset_names = [];
for n = 1:n_file
    dataset_names = [dataset_names '_' char(labels(n))];
end
filename = ['Saccades' dataset_names '.mat'];
save(fullfile(savedir, filename), 'stats', '-v7.3')
end