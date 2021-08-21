function [] = SOS_frf_free_fixed_all()
%% SOS_frf_free_fixed_all:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select head-free data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select head-fixed data', root, 'MultiSelect','off');
[FILE{3},PATH{3}] = uigetfile({'*.mat'}, 'Select body-fixed data', root, 'MultiSelect','off');

datasets = ["HeadFree", "HeadFixed", "BodyFixed"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'FRF_data','FUNC','U','N');
end
n_cond = ALL.HeadFree.N{1,3};

%% Virtual gaze (head-fixed body + body-fixed head)
clc
clearvars -except ALL FILE PATH n_cond

ALL.VirtualGaze = [];
ALL.VirtualGaze.FRF_data.IOFv = ALL.HeadFree.FRF_data.IOFv;
body_class = 'HeadFree';
for v = 1:n_cond
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).complex = ...
        ALL.(body_class).FRF_data.ref2body.grand_mean(v).complex + ...
        ALL.BodyFixed.FRF_data.ref2head.grand_mean(v).complex;
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).complex = sqrt( ...
        ALL.(body_class).FRF_data.ref2body.grand_std(v).complex.^(2) + ...
        ALL.BodyFixed.FRF_data.ref2head.grand_std(v).complex.^(2));
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).gain = ...
        abs(ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).complex);

    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).gain = sqrt( ...
        ALL.(body_class).FRF_data.ref2body.grand_std(v).gain.^(2) + ...
        ALL.BodyFixed.FRF_data.ref2head.grand_std(v).gain.^(2));
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).phase = ...
        rad2deg(angle(ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).complex));
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).phase = sqrt( ...
        ALL.(body_class).FRF_data.ref2body.grand_std(v).phase.^(2) + ...
        ALL.BodyFixed.FRF_data.ref2head.grand_std(v).phase.^(2));
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).error = ...
        abs((1 + 0i) - ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).complex);
    
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).error = sqrt( ...
        ALL.(body_class).FRF_data.ref2body.grand_std(v).error.^(2) + ...
        ALL.BodyFixed.FRF_data.ref2head.grand_std(v).error.^(2));
    
    % NaN
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).IO_coherence = nan*...
        ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).gain;
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).IO_coherence = nan*...
        ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).gain;
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).mag = nan*...
        ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).gain;
    ALL.VirtualGaze.FRF_data.ref2gaze.grand_std(v).mag = nan*...
        ALL.VirtualGaze.FRF_data.ref2gaze.grand_mean(v).gain;
end

%% Compare head free, head-fixed, body-fixed
set_names = ["HeadFree", "HeadFixed", "HeadFree", "HeadFree", "BodyFixed", "VirtualGaze"];
trf_names = ["ref2body", "ref2body", "ref2gaze", "ref2head", "ref2head", "ref2gaze"];
plot_names = ["mag", "gain", "phase", "error", "IO_coherence"];
yLines = [nan nan 0 1 nan];
n_set = length(trf_names);
n_plot = length(plot_names);
cond = ALL.HeadFree.U{1,3}{1};
cc = [0.9 0 0 ; 1 0.6 0.1 ; 0.5 0.3 1 ; 0 0.4 1 ; 0 0.8 0.2 ; 0.2 0.8 1];
% cc = distinguishable_colors(n_set);

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

set(ax(1,1:end),'YLim',[0 65])
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

%% Compare head free, head-fixed, body-fixed speeds
clc
clearvars -except ALL FILE PATH

set_names = ["HeadFree", "HeadFixed", "HeadFree", "HeadFree", "BodyFixed", "HeadFree"];
trf_names = ["ref2body", "ref2body", "ref2gaze", "ref2head", "ref2head", "body2head"];
plot_names = ["mag", "gain", "phase", "error", "IO_coherence"];
yLines = [nan nan 0 1 nan];
n_set = length(trf_names);
n_plot = length(plot_names);
n_cond = ALL.HeadFree.N{1,3};
cond = ALL.HeadFree.U{1,3}{1};
% cc = [0.9 0 0 ; 1 0.6 0.1; 0.5 0.3 1 ; 0 0.4 1 ; 0 0.8 0.2];
cc = flipud(distinguishable_colors(n_cond+6));

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 n_set*2 n_plot*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_plot,n_set);
for n = 1:n_set
    for p = 1:n_plot
        subI = n + (0:n_plot-1)*n_set;
        ax(p,n) = subplot(n_plot,n_set,subI(p)); hold on
        if p == 1
            title(trf_names(n) + "_" + set_names(n), 'interpreter', 'none')
        end
        if ~isnan(yLines(p))
            yline(yLines(p), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        end
        for v = 1:n_cond
            [h.patch(p,v,n),h.line(p,v,n)] = PlotPatch(...
                    ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_mean(v).(plot_names(p)),...
                  	ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_std(v).(plot_names(p)), ...
                    ALL.(set_names(n)).FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
        end
    end
end
leg_label = string(cond);
leg = legend(squeeze(h.line(1,:,1)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.005 0.91 0.07 0.08];
leg.Title.String = 'Speed (Â°/s)';

linkaxes(ax, 'x')
for a = 3:size(ax,1)
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

set(ax(1,1:end),'YLim',[0 65])
set(ax(2,1:end-1),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-250 150])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end-1), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end-1), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)
% delete(h.patch)

%% Error transforms
clc
clearvars -except ALL FILE PATH

set_names = ["HeadFree", "HeadFree"];
trf_names = ["err2body", "err2head"];
plot_names = ["gain", "phase", "IO_coherence"];
yLines = [nan 0 nan];
n_set = length(trf_names);
n_plot = length(plot_names);
n_cond = ALL.HeadFree.N{1,3};
cond = ALL.HeadFree.U{1,3}{1};
% cc = [0.9 0 0 ; 1 0.6 0.1; 0.5 0.3 1 ; 0 0.4 1 ; 0 0.8 0.2];
cc = flipud(distinguishable_colors(n_cond+6));

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 n_set*2 n_plot*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_plot,n_set);
for n = 1:n_set
    for p = 1:n_plot
        subI = n + (0:n_plot-1)*n_set;
        ax(p,n) = subplot(n_plot,n_set,subI(p)); hold on
        if p == 1
            title(trf_names(n) + "_" + set_names(n), 'interpreter', 'none')
        end
        if ~isnan(yLines(p))
            yline(yLines(p), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        end
        for v = 1:n_cond
            [h.patch(p,v,n),h.line(p,v,n)] = PlotPatch(...
                    ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_mean(v).(plot_names(p)),...
                  	ALL.(set_names(n)).FRF_data.(trf_names(n)).grand_std(v).(plot_names(p)), ...
                    ALL.(set_names(n)).FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
        end
    end
end
leg_label = string(cond);
leg = legend(squeeze(h.line(1,:,1)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.005 0.91 0.07 0.08];
leg.Title.String = 'Speed (°/s)';

linkaxes(ax, 'x')
for a = 2:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(2,1:end),'YLim',[0 1])
set(ax(2,1:end),'YLim',[-250 150])
set(ax(3,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end-1), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end-1), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)
% delete(h.patch)

%% Time constant
clc
clearvars -except ALL FILE PATH

set_names = ["HeadFixed", "HeadFree", "HeadFree"];
trf_names = ["ref2body", "ref2body", "ref2gaze"];
cc = [1 0.6 0.1; 0.9 0 0 ; 0.5 0.3 1];
full_names = trf_names + "_" + set_names;
n_set = length(trf_names);
n_cond = ALL.HeadFree.N{1,3};
cond = ALL.HeadFree.U{1,3}{1};

for v = 1:n_cond
    for n = 1:n_set
        time_const(v).(full_names(n)) = 1000*ALL.(set_names(n)).FRF_data.(trf_names(n)).fly(v).time_constant';
        r2(v).(full_names(n)) = ALL.(set_names(n)).FRF_data.(trf_names(n)).fly(v).time_constant_r2';
        G(v).(full_names(n)) = n * ones(size(r2(v).(full_names(n))));
    end
end

time_const_temp = squeeze(struct2cell(time_const));
r2_temp = squeeze(struct2cell(r2));
G_temp = squeeze(struct2cell(G));
for v = 1:n_cond
    time_const_group(:,v) = cat(1, time_const_temp{:,v});
    r2_group(:,v) = cat(1, r2_temp{:,v});
    G_group(:,v) = cat(1, G_temp{:,v});
end
G = G_group;
  
fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 2*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,n_cond);
for v = 1:n_cond
    ax(1,v) = subplot(2,n_cond,v); hold on
        b = boxchart(time_const_group(:,v),'GroupByColor', G(:,v), 'MarkerStyle', '.');
        
    ax(2,v) = subplot(2,n_cond,v + n_cond); hold on
        b = boxchart(r2_group(:,v),'GroupByColor', G(:,v), 'MarkerStyle', '.');
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

%% All time constant
fig = figure (5) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5 2*2])
movegui(fig, 'center')
clear ax h

time_const_keep = time_const_group(:);
r2_keep = r2_group(:);
G_keep = G(:);
r2_check = r2_keep < 0.5;
disp(['Removing ' num2str(sum(r2_check,'all')) ' of ' num2str(numel(r2_check)) ' flies'])

time_const_keep = time_const_keep(~r2_check);
r2_keep = r2_keep(~r2_check);
G_keep = G_keep(~r2_check);

rng(1)
spread = 0.3;
jitter = rand(size(G_keep)) * spread - (spread/2);
% cc = jet(n_set);
ax(1,1) = subplot(2,1,1); hold on
    %b = boxchart(time_const_keep,'GroupByColor', G_keep, 'MarkerStyle', '.');
    bx = boxplot(time_const_keep, G_keep, ...
        'Width', 0.5, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 0.5);

    h = get(bx(5,:),{'XData','YData'});
    for c = 1:size(h,1)
       patch(h{c,1},h{c,2}, cc(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'k','LineWidth', 1);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    
    plot(G_keep + jitter, time_const_keep, '.', 'Color', [0.5 0.5 0.5 0.2], 'MarkerSize', 2)

ax(2,1) = subplot(2,1,2); hold on
    %b = boxchart(r2_const_keep,'GroupByColor', G_keep, 'MarkerStyle', '.');
    bx = boxplot(r2_keep, G_keep, ...
        'Width', 0.5, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 0.5);

    h = get(bx(5,:),{'XData','YData'});
    for c = 1:size(h,1)
       patch(h{c,1},h{c,2}, cc(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end

    set(findobj(ax(2),'tag','Median'), 'Color', 'k','LineWidth', 1);
    set(findobj(ax(2),'tag','Box'), 'Color', 'none');
    set(findobj(ax(2),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(2),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(2).Children = ax(2).Children([end 1:end-1]);
    
    plot(G_keep + jitter, r2_keep, '.', 'Color', [0.5 0.5 0.5 0.2], 'MarkerSize', 2)
    
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

%% Stats
[p,~,stats] = anova1(time_const_keep, G_keep);
% [p,tb,stats] = kruskalwallis(time_const_keep, G_keep);
[c,m] = multcompare(stats, 'alpha', 0.001);

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