function [] = SOS_frf_free_fixed_error()
%% SOS_frf_free_fixed_error:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL');

%% Error transforms
clc
clearvars -except ALL FILE PATH

set_names = ["HeadFree", "BodyFixed"];
trf_names = [ "err2head", "err2head"];
plot_names = ["mag", "gain", "phase", "error", "IO_coherence"];
yLines = [nan nan 0 1 nan];
n_set = length(trf_names);
n_plot = length(plot_names);
cond = ALL.HeadFree.U{1,3}{1};
n_cond = ALL.HeadFree.N{1,3};
cc = [0 0.4 1; 0 0.8 0.2 ];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 n_cond*2 n_plot*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_plot,n_cond);
for v = 1:n_cond
    subI = v + (0:n_plot-1)*n_cond;
    for p = 1:n_plot
        ax(p,v) = subplot(n_plot,n_cond,subI(p)); hold on
        if p == 1
            title([num2str(cond(v)) '°/s'], 'interpreter', 'none')
        end
        for n = 1:n_set
            if ~isnan(yLines(p))
                yline(yLines(p), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            end
            [h.patch(p,v,n),h.line(p,v,n)] = PlotPatch(...
                    ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_mean(v).(plot_names(p)),...
                  	ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_std(v).(plot_names(p)), ...
                    ALL.(set_names(n)).FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
        end
    end
end
leg_label = trf_names + "_" + set_names;
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.16 0.38 0.63 0.1];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°/s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 100])
set(ax(2,1:end),'YLim',[0 1.3])
set(ax(3,1:end),'YLim',[-250 150])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)
% delete(h.patch)
end