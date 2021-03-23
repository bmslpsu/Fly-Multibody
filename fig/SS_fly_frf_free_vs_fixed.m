function [] = SS_fly_frf_free_vs_fixed()
%% SS_fly_frf_free_vs_fixed:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[HeadFree_file,HeadFree_path] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

% [HeadFixed_file,HeadFixed_path] = uigetfile({'*.mat'}, ...
%     'Select head fixed data', root, 'MultiSelect','off');

[BodyFixed_file,BodyFixed_path] = uigetfile({'*.mat'}, ...
    'Select body fixed data', root, 'MultiSelect','off');

HeadFree = load(fullfile(HeadFree_path,HeadFree_file),'FRF_data','FLY_mean','FUNC','U','N');
% HeadFixed = load(fullfile(HeadFixed_path,HeadFixed_file),'FRF_data','FLY_mean','FUNC','U','N');
BodyFixed = load(fullfile(BodyFixed_path,BodyFixed_file),'FRF_data','FLY_mean','FUNC','U','N');

%% Compare head
clc
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file

trf_names = ["ref2gaze", "ref2head", "ref2head"];
cc = [0.5 0.3 1;0 0.8 0.2; 0 0.4 1];
n_cond = 1;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 5*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_cond);
for v = 1
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(5,n_cond,subI(1)); hold on
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOMag,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOMag, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(1,v,2),h.line(1,v,2)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(2)).IOMag,...
                  BodyFixed.FRF_data.grand_std.(trf_names(2)).IOMag, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(1,v,3),h.line(1,v,3)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(3)).IOMag,...
                  HeadFree.FRF_data.grand_std.(trf_names(3)).IOMag, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(2,v) = subplot(5,n_cond,subI(2)); hold on
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOGain,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOGain, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(2,v,2),h.line(2,v,2)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(2)).IOGain,...
                  BodyFixed.FRF_data.grand_std.(trf_names(2)).IOGain, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(2,v,3),h.line(2,v,3)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(3)).IOGain,...
                  HeadFree.FRF_data.grand_std.(trf_names(3)).IOGain, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(3,v) = subplot(5,n_cond,subI(3)); hold on
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOPhaseDiff,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOPhaseDiff, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(3,v,2),h.line(3,v,2)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(2)).IOPhaseDiff,...
                  BodyFixed.FRF_data.grand_std.(trf_names(2)).IOPhaseDiff, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(3,v,3),h.line(3,v,3)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(3)).IOPhaseDiff,...
                  HeadFree.FRF_data.grand_std.(trf_names(3)).IOPhaseDiff, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(4,v) = subplot(5,n_cond,subI(4)); hold on
        yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(4,v,1),h.line(4,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOFRF_error,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOFRF_error, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(4,v,2),h.line(4,v,2)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(2)).IOFRF_error,...
                  BodyFixed.FRF_data.grand_std.(trf_names(2)).IOFRF_error, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(4,v,3),h.line(4,v,3)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(3)).IOFRF_error,...
                  HeadFree.FRF_data.grand_std.(trf_names(3)).IOFRF_error, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(5,v) = subplot(5,n_cond,subI(5)); hold on
        [h.patch(5,v,1),h.line(5,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOCohr,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOCohr, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(5,v,2),h.line(5,v,2)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(2)).IOCohr,...
                  BodyFixed.FRF_data.grand_std.(trf_names(2)).IOCohr, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(5,v,3),h.line(5,v,3)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(3)).IOCohr,...
                  HeadFree.FRF_data.grand_std.(trf_names(3)).IOCohr, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
end
leg_label = trf_names;
leg_label(1) = leg_label(1) + "_head-free";
leg_label(2) = leg_label(2) + "_body-fixed";
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.07 0.95 0.83 0.04];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

% delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.5 15],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

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

set(ax(1,1:end),'YLim',[0 50])
set(ax(2,1:end),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-250 150])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)

%% Compare wings
clc
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file

trf_names = ["ref2body", "ref2wing", "ref2wing"];
% cc = [0.9 0 0 ; 0.9 0.1 0.9; 0.1 0.9 0.2];
cc = [0.9 0 0 ; 0.1 0.7 1; 0.2 0.9 0.7];
n_cond = 1;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 5*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_cond);
for v = 1
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(5,n_cond,subI(1)); hold on
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOMag,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOMag, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(1,v,2),h.line(1,v,2)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(2)).IOMag,...
                  HeadFree.FRF_data.grand_std.(trf_names(2)).IOMag, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(1,v,3),h.line(1,v,3)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(3)).IOMag,...
                  BodyFixed.FRF_data.grand_std.(trf_names(3)).IOMag, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(2,v) = subplot(5,n_cond,subI(2)); hold on
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOGain,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOGain, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(2,v,2),h.line(2,v,2)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(2)).IOGain,...
                  HeadFree.FRF_data.grand_std.(trf_names(2)).IOGain, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(2,v,3),h.line(2,v,3)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(3)).IOGain,...
                  BodyFixed.FRF_data.grand_std.(trf_names(3)).IOGain, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(3,v) = subplot(5,n_cond,subI(3)); hold on
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOPhaseDiff,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOPhaseDiff, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(3,v,2),h.line(3,v,2)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(2)).IOPhaseDiff,...
                  HeadFree.FRF_data.grand_std.(trf_names(2)).IOPhaseDiff, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(3,v,3),h.line(3,v,3)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(3)).IOPhaseDiff,...
                  BodyFixed.FRF_data.grand_std.(trf_names(3)).IOPhaseDiff, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(4,v) = subplot(5,n_cond,subI(4)); hold on
        yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(4,v,1),h.line(4,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOFRF_error,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOFRF_error, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(4,v,2),h.line(4,v,2)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(2)).IOFRF_error,...
                  HeadFree.FRF_data.grand_std.(trf_names(2)).IOFRF_error, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(4,v,3),h.line(4,v,3)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(3)).IOFRF_error,...
                  BodyFixed.FRF_data.grand_std.(trf_names(3)).IOFRF_error, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
              
    ax(5,v) = subplot(5,n_cond,subI(5)); hold on
        [h.patch(5,v,1),h.line(5,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOCohr,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOCohr, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
        [h.patch(5,v,2),h.line(5,v,2)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(2)).IOCohr,...
                  HeadFree.FRF_data.grand_std.(trf_names(2)).IOCohr, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(2,:), 0.7*cc(2,:), 0.2, 1);
        [h.patch(5,v,3),h.line(5,v,3)] = PlotPatch(BodyFixed.FRF_data.grand_mean.(trf_names(3)).IOCohr,...
                  BodyFixed.FRF_data.grand_std.(trf_names(3)).IOCohr, BodyFixed.FRF_data.IOFv, ...
                  1, 1, cc(3,:), 0.7*cc(3,:), 0.2, 1);
end
leg_label = trf_names;
leg_label(1) = leg_label(1) + "_head-free";
leg_label(2) = leg_label(2) + "_head-free";
leg_label(3) = leg_label(3) + "_body-fixed";
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.07 0.95 0.83 0.04];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

% delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.5 15],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

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

set(ax(1,1:end),'YLim',[0 50])
set(ax(2,1:end),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-350 150])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)

%% wing2body
clc
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file

trf_names = ["wing2body"];
cc = [0 0.2 0.7];
cc = [0 0 0];
n_cond = 1;

fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 3*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_cond);
for v = 1
    subI = v + (0:2)*n_cond;
    ax(1,v) = subplot(3,n_cond,subI(1)); hold on
        h.fly(:,1) = plot(HeadFree.FRF_data.IOFv, HeadFree.FRF_data.fly_all.(trf_names(1)).IOGain, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5);
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOGain,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOGain, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
              
    ax(2,v) = subplot(3,n_cond,subI(2)); hold on
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        h.fly(:,2) = plot(HeadFree.FRF_data.IOFv, HeadFree.FRF_data.fly_all.(trf_names(1)).IOPhaseDiff - 180, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5);
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOPhaseDiff - 180,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOPhaseDiff, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
              
    ax(3,v) = subplot(3,n_cond,subI(3)); hold on
        h.fly(:,3) = plot(HeadFree.FRF_data.IOFv, HeadFree.FRF_data.fly_all.(trf_names(1)).IOCohr, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5);
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(HeadFree.FRF_data.grand_mean.(trf_names(1)).IOCohr,...
                  HeadFree.FRF_data.grand_std.(trf_names(1)).IOCohr, HeadFree.FRF_data.IOFv, ...
                  1, 1, cc(1,:), 0.7*cc(1,:), 0.2, 1);
  
end
leg_label = trf_names;
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.07 0.95 0.83 0.04];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

delete(h.patch)
set(h.fly, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 7)
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.5 15],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 400])
set(ax(2,1:end),'YLim',[-100 220])
set(ax(3,1:end),'YLim',[0 1.1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)

end