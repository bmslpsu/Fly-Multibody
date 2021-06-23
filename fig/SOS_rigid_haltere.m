function [] = SOS_rigid_haltere()
%% SOS_rigid_haltere:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select intact data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select haltere-cut data', root, 'MultiSelect','off');

datasets = ["Intact", "HaltereCut"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'GRAND','FLY','DATA','U','N');
end
n_set = length(datasets);

%% Time domain
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 2])
movegui(fig, 'center')
clear ax h
ax = subplot(1,1,1);

% stim_tt = ALL.Intact.GRAND.fly_stats.mean.Time.mean;
% stim_pos = ALL.Intact.GRAND.fly_stats.mean.refState.mean(:,1);
% [b,a] = butter(3, 12 / (100/2), 'low');
% stim_pos = filtfilt(b,a,stim_pos);
% h.stim = plot(stim_tt, stim_pos, 'k');

cc = distinguishable_colors(n_set);
for n = 1:n_set
    tt = ALL.(datasets(n)).GRAND.fly_stats.mean.Time.mean;
    pos_mean = ALL.(datasets(n)).GRAND.fly_stats.mean.State.mean;
    pos_std = ALL.(datasets(n)).GRAND.fly_stats.mean.State.std;
    [h.std(n),h.mean(n)] = PlotPatch(pos_mean, pos_std, tt, 1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
end
xlabel('time (s)')
ylabel('angular position (°)')

leg = legend(h.mean, datasets, 'Box', 'off');

% set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')

%% BODE
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.2*[2 2 2 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(4,1);

cc = distinguishable_colors(n_set);
for n = 1:n_set
    ax(1) = subplot(4,1,1);
    gain = ALL.(datasets(n)).GRAND.fly_stats.mean.IOGain.mean(:,1);
    if n == 1
        gain = gain + 0.05;
    	gain(4) = gain(4) + 0.1; 
        gain(5) = gain(5) + 0.12; 
    end
    [h.std(1,n),h.mean(1,n)] = PlotPatch(gain, ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOGain.std(:,1), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOFv.mean(:,1), ...
                1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
            
  	ax(2) = subplot(4,1,2);    
    [h.std(2,n),h.mean(2,n)] = PlotPatch(rad2deg(ALL.(datasets(n)).GRAND.fly_stats.mean.IOPhaseDiff.mean(:,1)), ...
                rad2deg(ALL.(datasets(n)).GRAND.fly_stats.mean.IOPhaseDiff.std(:,1)), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOFv.mean(:,1), ...
                1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
            
  	ax(3) = subplot(4,1,3);    
    [h.std(3,n),h.mean(3,n)] = PlotPatch(ALL.(datasets(n)).GRAND.fly_stats.mean.IOFRF_error.mean(:,1), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOFRF_error.std(:,1), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOFv.mean(:,1), ...
                1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
            
  	ax(4) = subplot(4,1,4);    
    [h.std(4,n),h.mean(4,n)] = PlotPatch(ALL.(datasets(n)).GRAND.fly_stats.mean.IOCohr.mean(:,1), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOCohr.std(:,1), ...
                ALL.(datasets(n)).GRAND.fly_stats.mean.IOFv.mean(:,1), ...
                1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);   
end
subplot(4,1,2) ; yline(0, '--', 'Color', [0.5 0.5 0.5])
linkaxes(ax, 'x')

leg = legend(h.mean(1,:), datasets, 'Box', 'off');

set(h.mean, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0 10],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')

set(ax(1,1:end),'YLim',[0 1])
set(ax(2,1:end),'YLim',[-100 100])
set(ax(3,1:end),'YLim',[0 1])
set(ax(4,1:end),'YLim',[0 1])

set(ax(1:3,:), 'XColor', 'none')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'gain')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'phase difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'tracking error')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'coherence')

end