function [] = SOS_time()
%% SOS_time:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','GRAND','FLY','D','I','U','N')

%% FFT of all wavelengths by fly
clearvars -except DATA ALL GRAND FLY D I U N root
clc

%% FRF
pI = [1 2 3 6];
T = ["body", "head", "gaze", "wing"];
n_plot = length(pI);
cc = hsv(n_plot);
Fv = DATA.reference{1}.Fv;

time = GRAND.all(1).Time(:,:,1);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_plot,1);
v = 1;
for n = 1:n_plot
    IOFv = GRAND.fly_stats(v).mean.IOFv.mean;
    subI = n + (0:2)*v;
    ax(n,1) = subplot(n_plot,1,n); hold on ; title(T(n))
        if (n==1) || (n==3)
            plot(time, GRAND.all_trial(v).refState.median(:,1), 'k', 'LineWidth', 1)
        end
        %plot(time, squeeze(GRAND.all(v).State(:,pI(n),:)), 'LineWidth', 0.5)
        [h.patch(1,n),h.line(1,n)] = PlotPatch(GRAND.all_trial(v).State.median(:,pI(n)),...
                  GRAND.all_trial(v).State.std(:,pI(n)), time, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

end
% set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
% set(ax, 'LineWidth', 1.5, 'FontSize', 12, 'XLim', [0.5 16],...
%     'XGrid','on','YGrid','on','box','on')
% 
% linkaxes(ax,'x')
% % linkaxes(ax(1,:),'y')
% linkaxes(ax(2,:),'y')
% linkaxes(ax(3,:),'y')
% 
% YLabelHC = get(ax(3,:), 'XLabel');
% set([YLabelHC{:}], 'String', 'Frequency (Hz)')
% 
% YLabelHC = get(ax(1,1), 'YLabel');
% set([YLabelHC], 'String', 'Gain (°/°)')
% YLabelHC = get(ax(2,1), 'YLabel');
% set([YLabelHC], 'String', 'Phase Difference (°)')
% YLabelHC = get(ax(3,1), 'YLabel');
% set([YLabelHC], 'String', 'Coherence')
% 
% set(ax(1,1:3),'YLim',[0 1.1])
% % set(ax(2), 'YLim', [0 0.5])
% % set(ax(4:6),'YLim',[-360 120])
% % set(ax(7:9),'YLim',[0 1.05])
% % set(ax(1:4),'XTickLabels',[])
% 
% % set(ax,'XScale','log')
% align_Ylabels(fig)


%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end