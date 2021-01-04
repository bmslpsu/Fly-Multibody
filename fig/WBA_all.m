function [] = WBA_all()
%% WBA_all: ncompare WBA of multiple data sets
root = 'E:\DATA\Magno_Data\Multibody';
[FILES,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');

n_file = length(FILES);
ALL = cell(n_file,1);
for n = 1:n_file
    ALL{n} = load(fullfile(PATH,FILES{n}),'DATA','D','I','U','N');
    filedata = textscan(char(FILES{n}), '%s', 'delimiter', '_');
    filedata = filedata{1};
    ALL{n}.name = [filedata{1} '_' filedata{2}];
end

%% WBA stats
clearvars -except ALL n_file FILES
clc

stats = cell(n_file,1);
for n = 1:n_file
    stats{n} = val_fly_stats(ALL{n});
end

%% Plot FFT by spatial wavelength
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 6])
ax = gobjects(3,N.wave);
pp = 1;
Fv = All.Fv{1}(:,1);
for w = [2:N.wave, 1]
    subI = pp + (0:2)*N.wave;
    ax(1,pp) = subplot(3,N.wave,subI(1)); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = All.wave_stats.head(w).mean;        
        plot(All.Fv{w}, All.head{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)
        
    pp = pp + 1;
end
set(ax, 'XLim', [0.5 40])
linkaxes(ax, 'x')
for p = 1:3
    linkaxes(ax(p,:), 'y')
end

%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
filename = ['Static_wave_fft_' clss '.mat'];
save(fullfile(savedir, filename), 'All','U','N','-v7.3')
end

%% Functions
function [stats] = val_fly_stats(all)
    DATA = all.DATA;
    [val_fly_group, val_group, ~] = findgroups(DATA{:,3}, DATA.fly);
    
    stats.lwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.lwing, ...
        'UniformOutput', false), val_fly_group);
    stats.rwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.rwing, ...
        'UniformOutput', false), val_fly_group);
    stats.dwba = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.dwba, ...
        'UniformOutput', false), val_fly_group);
    stats.head = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.head, ...
        'UniformOutput', false), val_fly_group);
    stats.body = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.body, ...
        'UniformOutput', false), val_fly_group);

    stats = structfun(@(x) splitapply(@(y) {y}, x, val_group), stats, 'UniformOutput', false);
    stats = structfun(@(x) cat(2,x{:}), stats, 'UniformOutput', false);

    stats.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
        stats, 'UniformOutput', false);

    fnames = string(fieldnames(stats.fly_stats));
    n_field = length(fnames);
    for f = 1:n_field
        for v = 1:all.N{1,3}
            stats.fly_mean.(fnames(f)){v} = cat(2, stats.fly_stats.(fnames(f))(:,v).mean);
        end
    end

    stats.val_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
        stats.fly_mean, 'UniformOutput', false);
end
