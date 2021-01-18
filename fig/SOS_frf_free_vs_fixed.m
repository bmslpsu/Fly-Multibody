function [] = SOS_frf_free_vs_fixed()
%% SOS_frf_free_vs_fixed:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[HeadFree_file,HeadFree_path] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

[HeadFixed_file,HeadFixed_path] = uigetfile({'*.mat'}, ...
    'Select head fixed data', root, 'MultiSelect','off');

[BodyFixed_file,BodyFixed_path] = uigetfile({'*.mat'}, ...
    'Select body fixed data', root, 'MultiSelect','off');

HeadFree = load(fullfile(HeadFree_path,HeadFree_file),'FRF_data','FUNC','U','N');
HeadFixed = load(fullfile(HeadFixed_path,HeadFixed_file),'FRF_data','FUNC','U','N');
BodyFixed = load(fullfile(BodyFixed_path,BodyFixed_file),'FRF_data','FUNC','U','N');

%%
clc
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file

trf_names = ["ref2body", "ref2body", "ref2gaze", "ref2head", "ref2head"];
n_cond = HeadFree.N{1,3};
n_plot = length(trf_names);
cc = hsv(n_plot);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_cond);
for v = 1:n_cond
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(5,n_cond,subI(1)); hold on
    title([num2str(HeadFree.U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).mag,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).mag, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(1,v,2),h.line(1,v,2)] = PlotPatch(HeadFixed.FRF_data.(trf_names(2)).grand_mean(v).mag,...
                  HeadFixed.FRF_data.(trf_names(2)).grand_std(v).mag, HeadFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(1,v,3),h.line(1,v,3)] = PlotPatch(HeadFree.FRF_data.(trf_names(3)).grand_mean(v).mag,...
                  HeadFree.FRF_data.(trf_names(3)).grand_std(v).mag, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        [h.patch(1,v,4),h.line(1,v,4)] = PlotPatch(HeadFree.FRF_data.(trf_names(4)).grand_mean(v).mag,...
                  HeadFree.FRF_data.(trf_names(4)).grand_std(v).mag, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(4,:), 0.7*cc(4,:), 0.2, 1);
        [h.patch(1,v,5),h.line(1,v,5)] = PlotPatch(BodyFixed.FRF_data.(trf_names(5)).grand_mean(v).mag,...
                  BodyFixed.FRF_data.(trf_names(5)).grand_std(v).mag, BodyFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(5,:), 0.7*cc(5,:), 0.2, 1);
              
    ax(2,v) = subplot(5,n_cond,subI(2)); hold on
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).gain,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).gain, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(2,v,2),h.line(2,v,2)] = PlotPatch(HeadFixed.FRF_data.(trf_names(2)).grand_mean(v).gain,...
                  HeadFixed.FRF_data.(trf_names(2)).grand_std(v).gain, HeadFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(2,v,3),h.line(2,v,3)] = PlotPatch(HeadFree.FRF_data.(trf_names(3)).grand_mean(v).gain,...
                  HeadFree.FRF_data.(trf_names(3)).grand_std(v).gain, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        [h.patch(2,v,4),h.line(2,v,4)] = PlotPatch(HeadFree.FRF_data.(trf_names(4)).grand_mean(v).gain,...
                  HeadFree.FRF_data.(trf_names(4)).grand_std(v).gain, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(4,:), 0.7*cc(4,:), 0.2, 1);
        [h.patch(2,v,5),h.line(2,v,5)] = PlotPatch(BodyFixed.FRF_data.(trf_names(5)).grand_mean(v).gain,...
                  BodyFixed.FRF_data.(trf_names(5)).grand_std(v).gain, BodyFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(5,:), 0.7*cc(5,:), 0.2, 1);
              
    ax(3,v) = subplot(5,n_cond,subI(3)); hold on
        yline(0, '--k', 'LineWidth', 1)
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).phase,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).phase, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(3,v,2),h.line(3,v,2)] = PlotPatch(HeadFixed.FRF_data.(trf_names(2)).grand_mean(v).phase,...
                  HeadFixed.FRF_data.(trf_names(2)).grand_std(v).phase, HeadFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(3,v,3),h.line(3,v,3)] = PlotPatch(HeadFree.FRF_data.(trf_names(3)).grand_mean(v).phase,...
                  HeadFree.FRF_data.(trf_names(3)).grand_std(v).phase, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        [h.patch(3,v,4),h.line(3,v,4)] = PlotPatch(HeadFree.FRF_data.(trf_names(4)).grand_mean(v).phase,...
                  HeadFree.FRF_data.(trf_names(4)).grand_std(v).phase, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(4,:), 0.7*cc(4,:), 0.2, 1);
        [h.patch(3,v,5),h.line(3,v,5)] = PlotPatch(BodyFixed.FRF_data.(trf_names(5)).grand_mean(v).phase,...
                  BodyFixed.FRF_data.(trf_names(5)).grand_std(v).phase, BodyFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(5,:), 0.7*cc(5,:), 0.2, 1);
              
    ax(4,v) = subplot(5,n_cond,subI(4)); hold on
        yline(1, '--k', 'LineWidth', 1)
        [h.patch(4,v,1),h.line(4,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).error,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).error, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(4,v,2),h.line(4,v,2)] = PlotPatch(HeadFixed.FRF_data.(trf_names(2)).grand_mean(v).error,...
                  HeadFixed.FRF_data.(trf_names(2)).grand_std(v).error, HeadFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(4,v,3),h.line(4,v,3)] = PlotPatch(HeadFree.FRF_data.(trf_names(3)).grand_mean(v).error,...
                  HeadFree.FRF_data.(trf_names(3)).grand_std(v).error, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        [h.patch(4,v,4),h.line(4,v,4)] = PlotPatch(HeadFree.FRF_data.(trf_names(4)).grand_mean(v).error,...
                  HeadFree.FRF_data.(trf_names(4)).grand_std(v).error, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(4,:), 0.7*cc(4,:), 0.2, 1);
        [h.patch(4,v,5),h.line(4,v,5)] = PlotPatch(BodyFixed.FRF_data.(trf_names(5)).grand_mean(v).error,...
                  BodyFixed.FRF_data.(trf_names(5)).grand_std(v).error, BodyFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(5,:), 0.7*cc(5,:), 0.2, 1);
              
    ax(5,v) = subplot(5,n_cond,subI(5)); hold on
        [h.patch(5,v,1),h.line(5,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).IO_coherence,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).IO_coherence, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(5,v,2),h.line(5,v,2)] = PlotPatch(HeadFixed.FRF_data.(trf_names(2)).grand_mean(v).IO_coherence,...
                  HeadFixed.FRF_data.(trf_names(2)).grand_std(v).IO_coherence, HeadFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(5,v,3),h.line(5,v,3)] = PlotPatch(HeadFree.FRF_data.(trf_names(3)).grand_mean(v).IO_coherence,...
                  HeadFree.FRF_data.(trf_names(3)).grand_std(v).IO_coherence, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
        [h.patch(5,v,4),h.line(5,v,4)] = PlotPatch(HeadFree.FRF_data.(trf_names(4)).grand_mean(v).IO_coherence,...
                  HeadFree.FRF_data.(trf_names(4)).grand_std(v).IO_coherence, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(4,:), 0.7*cc(4,:), 0.2, 1);
        [h.patch(5,v,5),h.line(5,v,5)] = PlotPatch(BodyFixed.FRF_data.(trf_names(5)).grand_mean(v).IO_coherence,...
                  BodyFixed.FRF_data.(trf_names(5)).grand_std(v).IO_coherence, BodyFixed.FRF_data.IOFv{v}, ...
                  1, 1, cc(5,:), 0.7*cc(5,:), 0.2, 1);
end
leg_label = trf_names;
leg_label(2) = leg_label(2) + "_head-fixed";
leg_label(5) = leg_label(5) + "_body-fixed";
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'horizontal');
leg.Position  = [0.17 0.95 0.63 0.05];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 11, 'LineWidth', 1.5)
set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid', 'on', 'YGrid', 'off', 'Box', 'on')
set(ax, 'XTick', [0.1, 1 10])
set(ax,'XScale','log')

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

set(ax(1,1:end),'YLim',[0 65])
set(ax(2,1:end-1),'YLim',[0 1.1])
set(ax(3,1:end),'YLim',[-250 250])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax,'XScale','log')
% align_Ylabels(fig)

%% Time constant
time_const = [];
r2 = [];
time_const_all = [];
r2_all = [];
for v = 1:n_cond
    time_const(v).body_head_free = 1000*HeadFree.FRF_data.(trf_names(1)).fly(v).time_constant';
    time_const(v).body_head_fixed = 1000*HeadFixed.FRF_data.(trf_names(1)).fly(v).time_constant';
    time_const(v).gaze_head_free = 1000*HeadFree.FRF_data.(trf_names(3)).fly(v).time_constant';

    r2(v).body_head_free = HeadFree.FRF_data.(trf_names(1)).fly(v).time_constant_r2';
    r2(v).body_head_fixed = HeadFixed.FRF_data.(trf_names(1)).fly(v).time_constant_r2';
    r2(v).gaze_head_free = HeadFree.FRF_data.(trf_names(3)).fly(v).time_constant_r2';

    time_const_all(:,v) = [time_const(v).body_head_free ; time_const(v).body_head_fixed ; ...
                                time_const(v).gaze_head_free ];
    r2_all(:,v) = [r2(v).body_head_free ; r2(v).body_head_fixed ; r2(v).gaze_head_free ];
end
G = [1*ones(size(time_const(v).body_head_free)) ; 2*ones(size(time_const(v).body_head_fixed)) ; ...
        3*ones(size(time_const(v).gaze_head_free))];

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 2*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,n_cond);
for v = 1:n_cond
    ax(1,v) = subplot(2,n_cond,v); hold on
        b = boxchart(time_const_all(:,v),'GroupByColor', G, 'MarkerStyle', '.');
        
    ax(2,v) = subplot(2,n_cond,v + n_cond); hold on
        b = boxchart(r2_all(:,v),'GroupByColor', G, 'MarkerStyle', '.');
end
fnames = fieldnames(time_const);
leg = legend(fnames, 'Box', 'off', 'interpreter', 'none', 'Orientation', 'horizontal');
leg.Position = [0.23 0.95 0.58 0.05];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end
set(ax(1,:), 'YLim', [-80 0])
set(ax(2,:), 'YLim', [0 1])

set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'Box', 'off', 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Time constant (ms)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'R^{2}')

%% Save time constant data
head_free_filedata = textscan(HeadFree_file, '%s', 'delimiter', '._');
head_fixed_filedata = textscan(HeadFixed_file, '%s', 'delimiter', '._');
body_fixed_filedata = textscan(BodyFixed_file, '%s', 'delimiter', '._');
head_free_name = [];
head_fixed_name = [];
body_fixed_name = [];
for n = 2:5
    head_free_name = [head_free_name '_' char(head_free_filedata{1}(n))];
    head_fixed_name = [head_fixed_name '_' char(head_fixed_filedata{1}(n))]; 
    body_fixed_name = [body_fixed_name '_' char(body_fixed_filedata{1}(n))]; 
end
fname = ['TimeConstant' head_free_name head_fixed_name body_fixed_name]; 

root = 'E:\DATA\Magno_Data\Multibody';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'time_const', 'time_const_all', 'r2', 'r2_all');
end