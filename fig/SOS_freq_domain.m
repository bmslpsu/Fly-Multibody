function [] = SOS_freq_domain()
%% SOS_freq_domain:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N')

%% Reference & Body
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

Fv = GRAND.all(1).Fv(:,:,1);
IOFv = GRAND.all(1).IOFv(:,:,1);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,1);
cc = [0.9 0 0 ; 0 0.5 0.9; 0 0.7 0.3];
for v = 1:N.vel
    IOFv = GRAND.fly_stats(v).median.IOFv.median(:,1);
    ax(v,1) = subplot(N.vel,1,v); hold on ; title([ num2str(U.vel{1}(v)) '°/s'])
        plot(FUNC{v}.All.Fv, FUNC{v}.All.mag.dX, 'k', 'LineWidth', 1)
        %plot(Fv, GRAND.all_trial(v).refMag.median(:,1), 'm', 'LineWidth', 1)  
        [h.patch(1,v),h.line(1,v)] = PlotPatch(GRAND.fly_stats(v).median.IOMag.median(:,1),...
                  GRAND.fly_stats(v).median.IOMag.std(:,1), IOFv, 0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(2,v),h.line(2,v)] = PlotPatch(GRAND.fly_stats(v).median.IOMag.median(:,2),...
                  GRAND.fly_stats(v).median.IOMag.std(:,2), IOFv, 0, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(3,v),h.line(3,v)] = PlotPatch(GRAND.fly_stats(v).median.IOMag.median(:,3),...
                  GRAND.fly_stats(v).median.IOMag.std(:,2), IOFv, 0, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        ax(v).YLim(1) = -4;
end
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 11, 'LineWidth', 1.5)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.2 15])
set(ax, 'XGrid', 'on', 'YGrid', 'off')
set(ax, 'XScale', 'log')

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', {'Velocity', 'magnitude (°/s)'})

YLabelHC = get(ax(end,1), 'XLabel');
set([YLabelHC], 'String', 'Frequency (Hz)')

align_Ylabels(fig)


end