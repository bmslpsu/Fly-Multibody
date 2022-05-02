function [] = SOS_frf_add_mass()
%% SOS_frf_free_fixed_all:
root = 'E:\EXPERIMENTS\MAGNO\Experiment_add_mass\data\processed';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'on');
FILE = sort(string(FILE))';

clear ALL
for n = 1:length(FILE)
    fileinfo = strsplit(FILE(n), '_');
    ALL.values(n) = str2double(fileinfo(3));
  	ALL.names(n) = join(fileinfo(2:3), '_');
    ALL.(ALL.names(n)) = load(fullfile(PATH,FILE(n)),'FRF_data','FUNC','U','N');
end

%% Compare head free, head-fixed, body-fixed speeds
clc
clearvars -except ALL FILE PATH

set_names = ["SOS_100", "SOS_200", "SOS_300", "SOS_500", "SOS_900"];
trf_names = ["ref2body", "ref2body", "ref2body", "ref2body", "ref2body"];
plot_names = ["gain", "phase", "error", "IO_coherence"];
yLines = [nan 0 1 nan];
n_set = length(trf_names);
n_plot = length(plot_names);
cc = (distinguishable_colors(n_set));
n_cond = 1;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 n_cond*2 n_plot*1.8])
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
leg.Position  = [0.2321    0.6994    0.5438    0.0829];

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end-1),'YLim',[0 1])
set(ax(2,1:end),'YLim',[-450 50])
set(ax(3,1:end),'YLim',[0 1.5])
set(ax(4,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end-1), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end-1), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)
% delete(h.patch)

%% Save FRF & time constant data
fname = 'FRF_combined';
for n = 1:length(FILE)
    filedata = textscan(FILE{n}, '%s', 'delimiter', '._');
    temp_name = [];
    for k = 2:5
        temp_name = [temp_name '_' char(filedata{1}(k))];
    end
    fname = [fname temp_name];
end
root = 'E:\DATA\Magno_Data\Multibody';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'ALL', 'time_const', 'time_const_keep', ...
    'r2', 'r2_keep', 'G', 'G_keep');
end