function [] = SOS_reference_time()
%% SOS_time:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','FUNC','GRAND','FLY','D','I','U','N')

%% Reference
clearvars -except DATA ALL GRAND FLY D I U N root
clc

time = GRAND.all(1).Time(:,:,1);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 6 4])
movegui(fig, 'center')
clear ax h
ax = gobjects(N{1,3},1);
sos_label = fieldnames(U);
sos_label = sos_label{3};
for v = 1:N{1,3}
    ax(v) = subplot(N{1,3},1,v); hold on ; title([sos_label ' = ' num2str(U{1,3}{1}(v))])
        plot(time, GRAND.all_trial(v).refState.median(:,1), 'k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1.5, 'Color', 'none', 'FontSize', 10, 'XLim', [-0.5 20])
linkaxes(ax,'xy')

XLabelHC = get(ax(3), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

YLabelHC = get(ax, 'YLabel');
set([YLabelHC{:}], 'String', 'Position (°)')

set(ax(1:2),'XTickLabels',[])

align_Ylabels(fig)

end