function [] = SOS_frf_add_mass()
%% SOS_frf_add_mass:
root = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_52_add_mass\data';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'DATA','FUNC','GRAND','U','N')

%%
clc
clearvars -except FILE DATA ALL GRAND FLY FUNC D I U N root

pI = [1 2];
T = ["ref2body", "err2body"];

% 3300
shift_I = {6:9, 6:9};
shift_I2 = {9, 9};
phase_lim = [-150 -150];
phase_lim2 = [-300 -300];

% % 1600
% shift_I = {8:9, 8:9};
% shift_I2 = {6:7, 6:77};
% phase_lim = [-200 -200];
% phase_lim2 = [-50 -50];

% % 900
% shift_I = {8:9, 8:9};
% shift_I2 = {7, 7};
% phase_lim = [-200 -200];
% phase_lim2 = [-50 -50];

% % 100
% shift_I = {8:9, 8:9};
% shift_I2 = {9, 9};
% phase_lim = [-50 -50];
% phase_lim2 = [-50 0];

fI = (1:9)';

FRF_data = [];
FRF_data.Fv = DATA.reference{1}.Fv;
n_plot = length(pI);
for v = 1:N{1,3}
    for n = 1:n_plot
        % Frequency components
        FRF_data.IOFv{v} = GRAND.fly_stats(v).mean.IOFv.mean(fI);
        n_freq = length(GRAND.fly_stats(v).mean.IOFv.mean);
        
        % FRF properties
        try
            mag_all = squeeze(GRAND.fly_all(v).mean.IOMag(fI,pI(n),:));
            mag_med = GRAND.fly_stats(v).mean.IOMag.mean(fI,pI(n));
            mag_std = GRAND.fly_stats(v).mean.IOMag.std(fI,pI(n));
        catch
            mag_all = squeeze(GRAND.fly_all(v).mean.refIOMag(fI,pI(1),:));
            mag_med = GRAND.fly_stats(v).mean.refIOMag.mean(fI,pI(1));
            mag_std = GRAND.fly_stats(v).mean.refIOMag.std(fI,pI(1));
        end

        gain_all = squeeze(GRAND.fly_all(v).mean.IOGain(fI,pI(n),:));
        gain_med = GRAND.fly_stats(v).mean.IOGain.mean(fI,pI(n));
        gain_std = GRAND.fly_stats(v).mean.IOGain.std(fI,pI(n));

        phase_all = rad2deg(squeeze(GRAND.fly_all(v).circ_mean.IOPhaseDiff(:,pI(n),:)));
        shift_all = any((1:n_freq)' == shift_I{n},2) & (phase_all > phase_lim(n));
        phase_all(shift_all) = phase_all(shift_all) - 360;
        
        shift_all = any((1:n_freq)' == shift_I2{n},2) & (phase_all > phase_lim2(n));
        phase_all(shift_all) = phase_all(shift_all) - 360;
        phase_all = phase_all(fI,:);

        phase_med = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_mean(:,pI(n)));
        phase_std = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_std(:,pI(n)));
        shift_all = any((1:n_freq)'==shift_I{n},2) & (phase_med > phase_lim(n));
        phase_med(shift_all) = phase_med(shift_all) - 360;
        shift_all = any((1:n_freq)'==shift_I2{n},2) & (phase_med > phase_lim2(n));
        phase_med(shift_all) = phase_med(shift_all) - 360;
        phase_med = phase_med(fI,:);
        phase_std = phase_std(fI,:);

        time_diff_all = 1000 * (phase_all ./360) .* (1 ./FRF_data.IOFv{v}(fI));
        time_diff_med = 1000 * (phase_med ./360) .* (1 ./FRF_data.IOFv{v}(fI));
        time_diff_std = 1000 * (phase_std ./360) .* (1 ./FRF_data.IOFv{v}(fI));

        cohr_all = squeeze(GRAND.fly_all(v).mean.Cohr(:,pI(n),:));
        cohr_med = GRAND.fly_stats(v).mean.Cohr.mean(:,pI(n));
        cohr_std = GRAND.fly_stats(v).mean.Cohr.std(:,pI(n));
        
      	IOcohr_all = squeeze(GRAND.fly_all(v).mean.IOCohr(fI,pI(n),:));
        IOcohr_med = GRAND.fly_stats(v).mean.IOCohr.mean(fI,pI(n));
        IOcohr_std = GRAND.fly_stats(v).mean.IOCohr.std(fI,pI(n));
        
        % Fly means
        FRF_data.(T(n)).fly(v).mag = mag_all;
        FRF_data.(T(n)).fly(v).gain = gain_all;
        FRF_data.(T(n)).fly(v).phase = phase_all;
        FRF_data.(T(n)).fly(v).time_diff = time_diff_all;
        FRF_data.(T(n)).fly(v).coherence = cohr_all;
        FRF_data.(T(n)).fly(v).IO_coherence = IOcohr_all;
        for f = 1:N.fly
            FRF_data.(T(n)).fly(v).complex(:,f) = gain_all(:,f) .* ...
                (cosd(phase_all(:,f)) + 1i*sind(phase_all(:,f)));
            FRF_data.(T(n)).fly(v).error(:,f) = abs((1 + 0*1i) - FRF_data.(T(n)).fly(v).complex(:,f));
            
            X = FRF_data.IOFv{v}(fI);
            Y = phase_all(:,f);
            [fitresult, gof] = time_constant_fit(X, Y, false);
            FRF_data.(T(n)).fly(v).time_constant(1,f) = fitresult.b;
            FRF_data.(T(n)).fly(v).time_constant_r2(1,f) = gof.rsquare;
        end
        
        % Grand means
        FRF_data.(T(n)).grand_mean(v).mag = mag_med;
        FRF_data.(T(n)).grand_mean(v).gain = gain_med;
        FRF_data.(T(n)).grand_mean(v).phase = phase_med;
        FRF_data.(T(n)).grand_mean(v).time_diff = time_diff_med;
        FRF_data.(T(n)).grand_mean(v).coherence = cohr_med;
        FRF_data.(T(n)).grand_mean(v).IO_coherence = IOcohr_med;
        FRF_data.(T(n)).grand_mean(v).complex = gain_med .* (cosd(phase_med) + 1i*sind(phase_med));
        FRF_data.(T(n)).grand_mean(v).error = mean(FRF_data.(T(n)).fly(v).error, 2);
        FRF_data.(T(n)).grand_mean(v).time_constant = mean(FRF_data.(T(n)).fly(v).time_constant, 2);
        FRF_data.(T(n)).grand_mean(v).time_constant_r2 = mean(FRF_data.(T(n)).fly(v).time_constant_r2, 2);
        
        % Grand STD's
        FRF_data.(T(n)).grand_std(v).mag = mag_std;
        FRF_data.(T(n)).grand_std(v).gain = gain_std;
        FRF_data.(T(n)).grand_std(v).phase = phase_std;
        FRF_data.(T(n)).grand_std(v).time_diff = time_diff_std;
        FRF_data.(T(n)).grand_std(v).coherence = cohr_std;
    	FRF_data.(T(n)).grand_std(v).IO_coherence = IOcohr_std;
        FRF_data.(T(n)).grand_std(v).error = std(FRF_data.(T(n)).fly(v).error, [], 2);
        FRF_data.(T(n)).grand_std(v).complex = std(FRF_data.(T(n)).fly(v).complex, [], 2);
        FRF_data.(T(n)).grand_std(v).time_constant = std(FRF_data.(T(n)).fly(v).time_constant, [], 2);
        FRF_data.(T(n)).grand_std(v).time_constant_r2 = std(FRF_data.(T(n)).fly(v).time_constant_r2, [], 2);
    end
end

%% FRF: one condition
cc = hsv(n_plot);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.2*n_plot 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(4,n_plot);
v = 1;
for n = 1:n_plot
    subI = n + (0:4)*n_plot;
    ax(1,n) = subplot(5,n_plot,subI(1)); hold on ; title(T(n), 'interpreter', 'none')
        %plot(FRF_data.IOFv{v}, GRAND.fly_stats(v).mean.refIOMag.mean(:,1), '*-', 'Color', 'k', 'LineWidth', 0.5)
        plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).mag, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(1,n),h.line(1,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).mag,...
                  FRF_data.(T(n)).grand_std(v).mag, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(2,n) = subplot(5,n_plot,subI(2)); hold on
        plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).gain, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(2,n),h.line(2,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).gain,...
                  FRF_data.(T(n)).grand_std(v).gain, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(3,n) = subplot(5,n_plot,subI(3)); hold on
        yline(0, '--k')
        plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).phase, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(3,n),h.line(3,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).phase,...
                  FRF_data.(T(n)).grand_std(v).phase, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(4,n) = subplot(5,n_plot,subI(4)); hold on
        yline(1, '--k')
        plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).error, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        %[h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).time_diff,...
                  %FRF_data.(T(n)).grand_std(v).time_diff, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
        [h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).error,...
                  FRF_data.(T(n)).grand_std(v).error, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(5,n) = subplot(5,n_plot,subI(5)); hold on
        plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).IO_coherence, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).IO_coherence,...
                  FRF_data.(T(n)).grand_std(v).IO_coherence, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    
        plot(FRF_data.Fv, FRF_data.(T(n)).fly(v).coherence, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
%         [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).coherence(2:end),...
%                   FRF_data.(T(n)).grand_std(v).coherence(2:end), FRF_data.Fv(2:end), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
end

set(h.line(1:4,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15')
set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 30],...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
linkaxes(ax(1,:),'y')
% linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
linkaxes(ax(4,:), 'y')
linkaxes(ax(5,:), 'y')

XLabelHC = get(ax(end,:), 'XLabel');
if n_plot == 1
    set([XLabelHC], 'String', 'Frequency (Hz)')
else
    set([XLabelHC{:}], 'String', 'Frequency (Hz)')
end

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (�)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (�/�)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (�)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(1,1:end),'YLim',[0 3.2])
% set(ax(2,1:end),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-600 100])
% set(ax(4,1:end),'YLim',300*[-1 1])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])

set(ax,'XScale','log')
align_Ylabels(fig)

%% FRF: all conditions
cc = hsv(N{1,3});

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.2*n_plot 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(4,n_plot);
fly_alpha = 0.3;
for v = 1:N{1,3}
    for n = 1:n_plot
        subI = n + (0:4)*n_plot;
        ax(1,n) = subplot(5,n_plot,subI(1)); hold on ; title(T(n), 'interpreter', 'none')
            %plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).mag, 'Color', [cc(v,:) fly_alpha], 'LineWidth', 0.5)
            [h.patch(1,n,v),h.line(1,n,v)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).mag,...
                      FRF_data.(T(n)).grand_std(v).mag, FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);

        ax(2,n) = subplot(5,n_plot,subI(2)); hold on
            %plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).gain, 'Color', [cc(v,:) fly_alpha], 'LineWidth', 0.5)
            [h.patch(2,n,v),h.line(2,n,v)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).gain,...
                      FRF_data.(T(n)).grand_std(v).gain, FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);

        ax(3,n) = subplot(5,n_plot,subI(3)); hold on
            yline(0, '--k')
            %plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).phase, 'Color', [cc(v,:) fly_alpha], 'LineWidth', 0.5)
            [h.patch(3,n,v),h.line(3,n,v)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).phase,...
                      FRF_data.(T(n)).grand_std(v).phase, FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);

        ax(4,n) = subplot(5,n_plot,subI(4)); hold on
            yline(1, '--k')
            %plot(FRF_data.IOFv{v}, FRF_data.(T(n)).fly(v).error, 'Color', [cc(v,:) fly_alpha], 'LineWidth', 0.5)
            [h.patch(4,n,v),h.line(4,n,v)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).error,...
                      FRF_data.(T(n)).grand_std(v).error, FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);

        ax(5,n) = subplot(5,n_plot,subI(5)); hold on
            %plot(FRF_data.Fv, FRF_data.(T(n)).fly(v).coherence, 'Color', [cc(v,:) fly_alpha], 'LineWidth', 0.5)
            [h.patch(5,n,v),h.line(5,n,v)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).coherence(2:end),...
                      FRF_data.(T(n)).grand_std(v).coherence(2:end), FRF_data.Fv(2:end), 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
                  
%             [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(T(n)).grand_mean(v).IO_coherence,...
%                   FRF_data.(T(n)).grand_std(v).IO_coherence, FRF_data.IOFv{v}, 1, 1, cc(v,:), 0.7*cc(v,:), 0.2, 1);
                  
    end
end
if N{1,3} > 1
    leg = legend(squeeze(h.line(5,end,:)), string(U{1,3}{1}), ...
        'Orientation', 'horizontal', 'Box', 'off');
    leg.Title.String = 'Stimulus speed (�/s)';
    leg.Position = [0.39 0.96 0.18 0.04];
end

% set(h.line(1:4,:,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 11')
% set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 20],...
%     'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
% set(ax, 'XTick', [0.1, 1 10])

set(h.line(1:4,:,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.2 25],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
linkaxes(ax(1,:),'y')
% linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
linkaxes(ax(4,:), 'y')
linkaxes(ax(5,:), 'y')

XLabelHC = get(ax(end,:), 'XLabel');
if n_plot == 1
    set([XLabelHC], 'String', 'Frequency (Hz)')
    set(ax(2,1:end),'YLim',[0 1.1])
else
    set([XLabelHC{:}], 'String', 'Frequency (Hz)')
    set(ax(2,1:end-2),'YLim',[0 1.1])
end

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (�/s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (�/�)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (�)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(1,1:end),'YLim',[0 80])
% set(ax(2,1:end-1),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-300 200])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
% set(ax(1:end-1,:), 'XTickLabel', [])
% set(ax(:,2:end-2), 'YTickLabels', [])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% align_Ylabels(fig)

%% Save FRF data
filedata = textscan(FILE, '%s', 'delimiter', '_');
dataset_name = [];
for n = 1:length(filedata{1})-1
    dataset_name = [dataset_name '_' char(filedata{1}(n))];
end
fname = ['FRF' dataset_name];
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'FRF_data', 'FUNC', 'U', 'N');

end