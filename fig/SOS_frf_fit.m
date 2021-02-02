function [] = SOS_frf_fit()
%% SOS_frf_fit:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[HeadFree_file,HeadFree_path] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');
HeadFree = load(fullfile(HeadFree_path,HeadFree_file),'FRF_data','FUNC','U','N');

%%
clearvars -except HeadFree
clc
v = 2;
Ts = 1/100;
IOFv = HeadFree.FRF_data.IOFv{v};
IOFv_rad = IOFv * 2*pi;
complex_response = HeadFree.FRF_data.ref2body.grand_mean(v).complex;

h = idfrd(complex_response, IOFv_rad, Ts);

plot(IOFv_rad, abs(complex_response), '-k*')

%% All speeds
trf_names = "head2body";
n_cond = HeadFree.N{1,3};
cc = colorcube(n_cond);
cc = [0 0.7 0.7 ; 0.9 0.9 0; 0.7 0 0.7];

all_phase_mean = mean(mean(cat(2,HeadFree.FRF_data.(trf_names(1)).grand_mean(:).phase),2));

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*1 3*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(3,1);
for v = 1:n_cond
    subI = 1:n_cond;
    gain = HeadFree.FRF_data.(trf_names(1)).grand_mean(v).gain;
    [~, h2b_switchI(v)] = min(abs(gain - 1));
    h2b_switch_gain(v) = gain(h2b_switchI(v));
    freq_switch(v) = HeadFree.FRF_data.IOFv{v}(h2b_switchI(v));
    ax(1,v) = subplot(3,1,subI(1)); hold on
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).gain,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).gain, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
              
    ax(2,v) = subplot(3,1,subI(2)); hold on
        %yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).phase,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).phase, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
              
    ax(3,v) = subplot(3,1,subI(3)); hold on
        [h.patch(3,v,1),h.line(3,v,1)] = PlotPatch(HeadFree.FRF_data.(trf_names(1)).grand_mean(v).IO_coherence,...
                  HeadFree.FRF_data.(trf_names(1)).grand_std(v).IO_coherence, HeadFree.FRF_data.IOFv{v}, ...
                  1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
end
leg = legend(h.line(1,:), string(HeadFree.U.vel{1}), 'Box', 'off');
linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.2 20],...
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

set(ax(1,1:end),'YLim',[0 40])
set(ax(2,1:end),'YLim',[-150 0])
set(ax(3,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(1:end-1,:), 'XColor', 'none')

set(ax,'XScale','log')
end