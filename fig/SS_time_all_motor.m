function [] = SS_time_all_motor()
%% SS_time_all_motor:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','Replay','GRAND','D','I','U','N')

%% Reference & Body
clearvars -except FUNC Replay DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 1.3*N.freq])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.freq,1);
cc = [0.9 0 0 ; 0 0.4 1];
for v = 1:N.freq
    ax(v,1) = subplot(N.freq,1,v); hold on ; %title([ num2str(U.freq{1}(v)) 'hz'])
        %plot(FUNC{v}.All.time, FUNC{v}.All.X, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
        %plot(Replay{1}.time, Replay{1}.pos.body_sine(:,v), 'k', 'LineWidth', 1)
        %plot(time, GRAND.all_trial(v).refState.median(:,1), 'm', 'LineWidth', 1)
        %scale = max(GRAND.fly_stats(v).mean.State.mean(:,1)) / max(GRAND.all_trial(v).refState.mean(:,1));
        scale = 1;
        plot(time, squeeze(GRAND.all(v).State(:,1,:)), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
%         [h.patch(1,v),h.line(1,v)] = PlotPatch(scale*GRAND.all_trial(v).refState.mean(:,1),...
%                   GRAND.all_trial(v).State.std(:,1), time, 0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);     
        [h.patch(2,v),h.line(2,v)] = PlotPatch(GRAND.fly_stats(v).mean.State.mean(:,1),...
                  GRAND.fly_stats(v).mean.State.std(:,1), time, 1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
%         set(ax(v), 'YLim', max(abs(ax(v).YLim))*[-1 1])
end
set(ax, 'Color', 'none', 'LineWidth', 1.5, 'FontSize', 10, 'XLim', [-0.2 10], 'XTick', 0:2:20)
set(ax(1:end-1,:), 'XColor', 'none')
set(ax, 'YLim', 10*[-1 1])

% set(ax, 'YLim', 10*[-1 1])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Angular position (°)')

YLabelHC = get(ax(end,1), 'XLabel');
set([YLabelHC], 'String', 'Time (s)')

align_Ylabels(fig)

end