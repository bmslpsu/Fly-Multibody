function [] = SOS_time_fly()
%% SOS_time:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N')

%% FLY MEAN: Reference, Body, Head, Gaze, WBA
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 5])
movegui(fig, 'center')
clear ax h

flyI = 9;
tI = [1 2 3 5];
yL = ["Body", "Head", "Gaze", "\DeltaWBA"] + " (°)";
cc = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1 ; 0 0.6 0.4];
n_plot = length(tI);
ax = gobjects(n_plot,N.vel);
for n = 1:n_plot
    for v = 1:N.vel
        subI = N.vel*(n-1) + v;
        ax(n,v) = subplot(n_plot,N.vel,subI); hold on
            if v == 1
                ylabel(yL(n))
            end

            if n == 1
                title([ num2str(U.vel{1}(v)) '°s^{-1}'])
                plot(FUNC{v}.All.time, FUNC{v}.All.X, 'k', 'LineWidth', 1)
            end

            [h.patch(n,v),h.line(n,v)] = PlotPatch(FLY.stats(flyI,v).State.mean(:,tI(n)),...
                      FLY.stats(flyI,v).State.std(:,tI(n)), ...
                      time, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
            set(ax(n,v), 'YLim', (1 + max(abs(ax(n,v).YLim)))*[-1 1])
    end
    linkaxes(ax(n,:), 'y')
end
linkaxes(ax, 'x')
linkaxes(ax([1,3],:), 'y')
set(ax, 'Color', 'none', 'LineWidth', 1.5, 'FontSize', 10, 'XLim', [-0.5 20], 'XTick', 0:5:20)

YLabelHC = get(ax(end,:), 'XLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
set(ax(1:end-1,:), 'XColor', 'none')

align_Ylabels(ax(:,1))

%% TRIAL: Reference, Body, Head, Gaze, WBA
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 5])
movegui(fig, 'center')
clear ax h

flyI = 9;
trialI = 10;
tI = [1 2 3 5];
yL = ["Body", "Head", "Gaze", "\DeltaWBA"] + " (°)";
cc = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1 ; 0 0.6 0.4];
n_plot = length(tI);
ax = gobjects(n_plot,N.vel);
for n = 1:n_plot
    for v = 1:N.vel
        subI = N.vel*(n-1) + v;
        ax(n,v) = subplot(n_plot,N.vel,subI); hold on
            if v == 1
                ylabel(yL(n))
            end

            if n == 1
                title([ num2str(U.vel{1}(v)) '°s^{-1}'])
                plot(FUNC{v}.All.time, FUNC{v}.All.X, 'k', 'LineWidth', 1)
            end

            trial = FLY.all(flyI,v).State(:,tI(n),trialI);
            [h.patch(n,v),h.line(n,v)] = PlotPatch(trial,...
                      0*trial, ...
                      time, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
            set(ax(n,v), 'YLim', (1 + max(abs(ax(n,v).YLim)))*[-1 1])
    end
    linkaxes(ax(n,:), 'y')
end
linkaxes(ax, 'x')
linkaxes(ax([1,3],:), 'y')
set(ax, 'Color', 'none', 'LineWidth', 1.5, 'FontSize', 10, 'XLim', [-0.5 20], 'XTick', 0:5:20)

YLabelHC = get(ax(end,:), 'XLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
set(ax(1:end-1,:), 'XColor', 'none')

align_Ylabels(ax(:,1))


end