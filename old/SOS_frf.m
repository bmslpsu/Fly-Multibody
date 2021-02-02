function [] = SOS_frf()
%% SOS_frf:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','FUNC','GRAND','FLY','D','I','U','N')

%% FRF
clearvars -except DATA ALL GRAND FLY D I U N root
clc
% pI = [1 2 3 6 8 9 11 12];
% T = ["ref2body", "ref2head", "ref2gaze", "ref2wing", "head2body", "head2wing", "wing2body", "left2right"];

% pI = [1 2 3 8];
% T = ["ref2body", "ref2head", "ref2gaze", "head2body"];

pI = [1];
T = ["ref2body_fixed"];

% pI = [1 2 3];
% T = ["ref2head", "ref2wing", "head2wing"];

n_plot = length(pI);
cc = hsv(n_plot);
Fv = DATA.reference{1}.Fv;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.2*n_plot 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(4,n_plot);

shift_I = {3:7, 5, 6, 3:7};
phase_lim = [0 nan nan 0];

% shift_I = {7:9, 9, 9, 9, 9, 1:9, 3:7, 7, 6, 3:7};
% phase_lim = [0 nan nan 0 0 0 nan 0 -100 0 nan 0];
fI = (1:7)';
for v = 3
    for n = 1:n_plot
        IOFv = GRAND.fly_stats(v).mean.IOFv.mean;
        n_freq = length(IOFv);
        subI = n + (0:4)*n_plot;

        ax(1,n) = subplot(5,n_plot,subI(1)); hold on ; title(T(n), 'interpreter', 'none')
            mag_all = squeeze(GRAND.fly_all(v).mean.IOMag(:,pI(n),:));
            plot(IOFv(fI), mag_all((fI),:), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            mag_med = GRAND.fly_stats(v).mean.IOMag.mean(fI,pI(n));
            mag_std = GRAND.fly_stats(v).mean.IOMag.std(fI,pI(n));
            stim_med = GRAND.fly_stats(v).mean.refIOMag.mean(fI,1);
            plot(IOFv(fI), stim_med, 'k', 'Marker', '.','MarkerFaceColor', 'none', ...
                                            'MarkerSize', 15')
            [h.patch(1,n),h.line(1,n)] = PlotPatch(mag_med,...
                      mag_std, IOFv(fI), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
        ax(2,n) = subplot(5,n_plot,subI(2)); hold on
            %gain_all = squeeze(GRAND.all(v).IOGain(:,pI(n),:));
            gain_all = squeeze(GRAND.fly_all(v).mean.IOGain(:,pI(n),:));
            plot(IOFv(fI), gain_all((fI),:), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            gain_med = GRAND.fly_stats(v).mean.IOGain.mean(fI,pI(n));
            gain_std = GRAND.fly_stats(v).mean.IOGain.std(fI,pI(n));
            [h.patch(2,n),h.line(2,n)] = PlotPatch(gain_med,...
                      gain_std, IOFv(fI), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    
        ax(3,n) = subplot(5,n_plot,subI(3)); hold on
            plot([0 20], [0 0], '--k')
            
            %phase_all = rad2deg(squeeze(GRAND.all(v).IOPhaseDiff(:,pI(n),:)));
            phase_all = rad2deg(squeeze(GRAND.fly_all(v).circ_mean.IOPhaseDiff(:,pI(n),:)));
            shift_all = any((1:n_freq)' == shift_I{n},2) & (phase_all > phase_lim(n));
            phase_all(shift_all) = phase_all(shift_all) - 360;  
            plot(IOFv(fI), phase_all((fI),:), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
     
            phase_med = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_mean(:,pI(n)));
            phase_std = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_std(:,pI(n)));
            shift_all = any((1:n_freq)'==shift_I{n},2) & (phase_med > phase_lim(n));
            phase_med(shift_all) = phase_med(shift_all) - 360;  
            
            [h.patch(3,n),h.line(3,n)] = PlotPatch(phase_med(fI),...
                      phase_std(fI), IOFv(fI), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
        ax(4,n) = subplot(5,n_plot,subI(4)); hold on
            timediff_all = 1000 * (phase_all ./360) .* (1 ./IOFv);
            plot(IOFv(fI), timediff_all((fI),:), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            
            timediff_med = 1000 * (phase_med ./360) .* (1 ./IOFv);
            timediff_std = 1000 * (phase_std ./360) .* (1 ./IOFv);
            [h.patch(4,n),h.line(4,n)] = PlotPatch(timediff_med(fI),...
                      timediff_std(fI), IOFv(fI), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

        ax(5,n) = subplot(5,n_plot,subI(5)); hold on
            %cohr_all = squeeze(GRAND.all(v).Cohr(:,pI(n),:));
            cohr_all = squeeze(GRAND.fly_all(v).mean.Cohr(:,pI(n),:));
            plot(Fv, cohr_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            cohr_med = GRAND.fly_stats(v).mean.Cohr.mean(:,pI(n));
            cohr_std = GRAND.fly_stats(v).mean.Cohr.std(:,pI(n));
            [h.patch(5,n),h.line(5,n)] = PlotPatch(cohr_med,...
                      cohr_std, Fv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    end
end
% set(ax,'XScale','log')

set(h.line(1:4,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15')
set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid','on','YGrid','on','box','on')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
linkaxes(ax(1,:),'y')
% linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
linkaxes(ax(4,:), 'y')
linkaxes(ax(5,:), 'y')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Time difference (ms)')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 100])
set(ax(2,1:end),'YLim',[0 1.3])
set(ax(3,1:end),'YLim',[-300 200])
set(ax(4,1:end),'YLim',300*[-1 1])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])

% set(ax(2), 'YLim', [0 0.5])
% set(ax(4:6),'YLim',[-360 120])
% set(ax(7:9),'YLim',[0 1.05])
% set(ax(1:4),'XTickLabels',[])

set(ax,'XScale','log')
align_Ylabels(fig)


%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end