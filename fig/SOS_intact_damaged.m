function [] = SOS_intact_damaged()
%% SOS_intact_damaged:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[Intact_file,Intact_path] = uigetfile({'*.mat'}, ...
    'Select intact fly data', root, 'MultiSelect','off');

[Damaged_file,Damaged_path] = uigetfile({'*.mat'}, ...
    'Select damaged fly data', root, 'MultiSelect','off');

DATA{1} = load(fullfile(Intact_path,Intact_file),'FRF_data','FLY_mean','FUNC','U','N');
DATA{2} = load(fullfile(Damaged_path,Damaged_file),'FRF_data','FLY_mean','FUNC','U','N');
n_file = length(DATA);

%% Compare head
clc
clearvars -except DATA n_file

trf_names = ["wing2body"];
cc = [0.5 0.3 1; 0.1 0.8 0.4];
n_cond = 1;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 3*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(3,n_cond);
for v = 1
    subI = v + (0:4)*n_cond;
    ax(1,v) = subplot(3,n_cond,subI(1)); hold on
        for n = 1: n_file
            [h.patch(1,v,n),h.line(1,v,n)] = PlotPatch(DATA{n}.FRF_data.grand_mean.(trf_names(1)).IOGain,...
                      DATA{n}.FRF_data.grand_std.(trf_names(1)).IOGain, DATA{n}.FRF_data.IOFv{v}, ...
                      1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
        end
              
    ax(2,v) = subplot(3,n_cond,subI(2)); hold on
        yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        for n = 1: n_file
            [h.patch(2,v,n),h.line(2,v,n)] = PlotPatch(DATA{n}.FRF_data.grand_mean.(trf_names(1)).IOPhaseDiff,...
                      DATA{n}.FRF_data.grand_std.(trf_names(1)).IOPhaseDiff, DATA{n}.FRF_data.IOFv{v}, ...
                      1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
        end
              
    ax(3,v) = subplot(3,n_cond,subI(3)); hold on
        for n = 1: n_file
            [h.patch(3,v,n),h.line(3,v,n)] = PlotPatch(DATA{n}.FRF_data.grand_mean.(trf_names(1)).IOCohr,...
                      DATA{n}.FRF_data.grand_std.(trf_names(1)).IOCohr, DATA{n}.FRF_data.IOFv{v}, ...
                      1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
        end
end
leg_label = trf_names;
leg_label(1) = trf_names(1) + "_intact";
leg_label(2) = trf_names(1) + "_damaged";
leg = legend(squeeze(h.line(1,1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.07 0.95 0.83 0.04];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

% delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 15],...
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

% set(ax(1,1:end),'YLim',[0 50])
set(ax(2,1:end),'YLim',[-300 300])
set(ax(3,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
align_Ylabels(fig)

%% Stats
cmp = ["IOGain", "IOPhaseDiff"];
P = [];
H = [];
alpha = 0.005;
for c = 1:length(cmp)
    Intact = DATA{1}.FRF_data.fly_all.wing2body.(cmp(c));
    Damaged = DATA{2}.FRF_data.fly_all.wing2body.(cmp(c));

    n_freq = size(Intact,1);
    P.(cmp(c)) = nan(n_freq,1);
    for f = 1:n_freq
        I = Intact(f,:);
        D = Damaged(f,:);
        I = I(~isnan(I));
        D = D(~isnan(D));

        T = [I' ; D'];
        G = [ones(length(I),1) ; 2*ones(length(D),1)];

        %[~,P.(cmp(c))(f,1)] = ttest(I',D');
        [P.(cmp(c))(f,1)] = anova1(T, G, 'off');
        H.(cmp(c))(f,1) = P.(cmp(c))(f,1) < alpha;
    end
end

end