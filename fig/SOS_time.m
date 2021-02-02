function [] = SOS_time()
%% SOS_time:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N')

%% Reference & Body
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,1);
cc = [0.9 0 0 ; 0 0.7 0.3];
for v = 1:N.vel
    ax(v,1) = subplot(N.vel,1,v); hold on ; title([ num2str(U.vel{1}(v)) '°/s'])
        plot(FUNC{v}.All.time, FUNC{v}.All.X, 'k', 'LineWidth', 1)
        %plot(time, GRAND.all_trial(v).refState.median(:,1), 'm', 'LineWidth', 1)
        %[h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.all_trial(v).State.median(:,pI(1)),...
                  %GRAND.all_trial(v).State.std(:,pI(1)), time, 0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);     
        [h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.fly_stats(v).median.State.mean(:,1),...
                  GRAND.fly_stats(v).median.State.std(:,1), time, 1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        set(ax(v), 'YLim', max(abs(ax(v).YLim))*[-1 1])
end
set(ax, 'Color', 'none', 'LineWidth', 1.5, 'FontSize', 10, 'XLim', [-0.2 20], 'XTick', 0:2:20)

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', 'Angular position (°)')

YLabelHC = get(ax(end,1), 'XLabel');
set([YLabelHC], 'String', 'Time (s)')

align_Ylabels(fig)

%% Head
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,1);
cc = [0 0.4 0.8];
for v = 1:N.vel
    ax(v,1) = subplot(N.vel,1,v); hold on ; title([ num2str(U.vel{1}(v)) '°/s'])
        %[h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.all_trial(v).State.median(:,2),...
                  %GRAND.all_trial(v).State.std(:,2), time, 0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.fly_stats(v).median.State.median(:,2),...
                  GRAND.fly_stats(v).median.State.std(:,2), time, 1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
    set(ax(v), 'YLim', (1 + max(abs(ax(v).YLim)))*[-1 1])
end
set(ax, 'Color', 'none', 'LineWidth', 1.5, 'FontSize', 10, 'XLim', [-0.2 20], 'XTick', 0:2:20)

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', 'Angular position (°)')

YLabelHC = get(ax(end,1), 'XLabel');
set([YLabelHC], 'String', 'Time (s)')

align_Ylabels(fig)

end