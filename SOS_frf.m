function [] = SOS_frf()
%% SOS_magnitude:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','GRAND','FLY','FUNC','D','I','U','N')

%% FRF
clearvars -except DATA ALL GRAND FLY FUNC D I U N root
clc

% lwI = 4;
% rwI = 5;
% dwI = 6;
% h2b = 8;
% h2w = 9;
% w2b = 11;

pI = [1 2 3 6 8 9 11 13];
T = ["ref2body", "ref2head", "ref2gaze", "ref2wing", ...
    "head2body", "wing2head", "wing2body","left2right"];
n_plot = length(pI);
cc = hsv(n_plot);
Fv = DATA.reference{1}.Fv;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
movegui(fig, 'center')
clear ax h
v = 1;
ax = gobjects(3,n_plot);
for n = 1:n_plot
    IOFv = GRAND.fly_stats(v).mean.IOFv.mean;
    subI = n + (0:2)*n_plot;
    ax(1,n) = subplot(3,n_plot,subI(1)); hold on ; title(T(n))
        grand_med = GRAND.fly_stats(v).mean.IOGain.mean(:,pI(n));
        grand_std = GRAND.fly_stats(v).std.IOGain.std(:,pI(n));
        [h.patch(1,n),h.line(1,n)] = PlotPatch(grand_med,...
                  grand_std, IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(2,n) = subplot(3,n_plot,subI(2)); hold on
        plot([0 10], [0 0], '--k')
        %phase = grand_mean_phase(:,pI(n));
        %phase(phase > phase_lim(n)) = phase(phase > phase_lim(n)) - 360;
        %phase(any((1:N.freq)'==shift_I{n},2)) = phase(any((1:N.freq)'==shift_I{n},2)) - 360;       
        grand_med = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_mean(:,pI(n)));
        grand_std = rad2deg(GRAND.fly_stats(v).circ_std.IOPhaseDiff.circ_mean(:,pI(n)));
        grand_med(grand_med > 20) = grand_med(grand_med > 20) - 360;
        [h.patch(2,n),h.line(2,n)] = PlotPatch(grand_med,...
                  grand_std, IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(3,n) = subplot(3,n_plot,subI(3)); hold on
        grand_med = GRAND.fly_stats(v).mean.Cohr.mean(:,pI(n));
        grand_std = GRAND.fly_stats(v).std.Cohr.mean(:,pI(n));
        [h.patch(3,n),h.line(3,n)] = PlotPatch(grand_med,...
                  grand_std, Fv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

end
set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontSize', 12, 'XLim', [1 16],...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')

linkaxes(ax,'x')
% linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
linkaxes(ax(3,:),'y')

YLabelHC = get(ax(3,:), 'XLabel');
set([YLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase Difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(1:3),'YLim',[0 1.05])
% set(ax(2), 'YLim', [0 0.5])
% set(ax(4:6),'YLim',[-360 120])
% set(ax(7:9),'YLim',[0 1.05])
% set(ax(1:4),'XTickLabels',[])

% set(ax,'XScale','log')

%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end