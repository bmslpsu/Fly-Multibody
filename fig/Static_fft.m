function [] = Static_fft()
%% Static_fft:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'DATA','D','I','U','N')

%% Average by fly
clearvars -except DATA N U
[wave_fly_group, wave_group, fly_group] = findgroups(DATA.wave, DATA.fly);

clss = 'position';
All.Fv = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.Fv, DATA.body, ...
    'UniformOutput', false), wave_fly_group);
All.body = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.mag.(clss), DATA.body, ...
    'UniformOutput', false), wave_fly_group);
All.head = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.mag.(clss), DATA.head, ...
    'UniformOutput', false), wave_fly_group);
All.dwba = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.mag.(clss), DATA.dwba, ...
    'UniformOutput', false), wave_fly_group);

All = structfun(@(x) splitapply(@(y) {y}, x, findgroups(wave_group)), All, 'UniformOutput', false);
All = structfun(@(x) cat(2,x{:}), All, 'UniformOutput', false);

All.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    All, 'UniformOutput', false);

fnames = string(fieldnames(All.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:N.wave
        All.fly_mean.(fnames(f)){v} = cat(2, All.fly_stats.(fnames(f))(:,v).mean);
        All.all_trial.(fnames(f)){v} = cat(2, All.(fnames(f)){:,v});
    end
end

All.wave_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    All.fly_mean, 'UniformOutput', false);

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
        plot(All.all_trial.Fv{w}, All.all_trial.head{w}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        %plot(All.fly_mean.Fv{w}, All.fly_mean.head{w}, 'b')
        plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)

        xline(12,'--r')
        xline(19,'--r')
        
    ax(2,pp) = subplot(3,N.wave,subI(2)); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = All.wave_stats.body(w).mean;        
        plot(All.all_trial.Fv{w}, All.all_trial.body{w}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)

        xline(12,'--r')
        xline(19,'--r')
        
    ax(3,pp) = subplot(3,N.wave,subI(3)); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = All.wave_stats.dwba(w).mean;        
        plot(All.all_trial.Fv{w}, All.all_trial.dwba{w}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)

        xline(12,'--r')
        xline(19,'--r')
        
    pp = pp + 1;
end
set(ax, 'XLim', [1 40])
linkaxes(ax, 'x')
for p = 1:3
    linkaxes(ax(p,:), 'y')
end

set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'Box', 'off')
set(ax(:,2:end), 'XColor', 'none', 'YColor', 'none')
set(ax(1:end-1,:), 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Head (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Body (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', '\DeltaWBA (°)')

XLabelHC = get(ax(3,1), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

set(ax, 'YLim', [-0.1 1])

%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
filename = ['Static_wave_fft_' clss '.mat'];
save(fullfile(savedir, filename), 'All','U','N','-v7.3')

end