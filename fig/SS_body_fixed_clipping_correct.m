function [] = SS_body_fixed_clipping_correct()
%% SS_body_fixed_clipping_correct:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select head-free data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select body-fixed data', root, 'MultiSelect','off');
datasets = ["HeadFree", "BodyFixed"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'GRAND','FUNC','U','N');
end

%% Compare head free, head-fixed, body-fixed
clc ; close all
clearvars -except ALL FILE PATH root

% set_names = ["BodyFixed"];
% state_names = ["head"];
% idx_state = [1];
% cc = [0 0.4 1];

set_names = ["HeadFree", "HeadFree", "BodyFixed"];
state_names = ["body", "head", "head"];
idx_state = [1 2 1];

n_set = length(state_names);
n_cond = ALL.HeadFree.N.freq;
freq = ALL.HeadFree.U.freq{1};

tt = ALL.(set_names(1)).GRAND.fly_stats(1).mean.Time.mean(:,idx_state(1));
fs = 1 / mean(diff(tt));
fit_range = (1*fs + 1):length(tt);
fft_frf = [];
FRF_fit = [];
for v = 1:n_cond
    for n = 1:n_set
        fft_frf.(set_names(n)).(state_names(n)).gain.mean(v) = ...
             	ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOGain.mean(:, idx_state(n));
        fft_frf.(set_names(n)).(state_names(n)).gain.std(v) = ...
                ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOGain.std(:, idx_state(n));

        fft_frf.(set_names(n)).(state_names(n)).phase.mean(v) = ...
                rad2deg(ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOPhaseDiff.mean(:, idx_state(n)));
        fft_frf.(set_names(n)).(state_names(n)).phase.std(v) = ...
                rad2deg(ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOPhaseDiff.std(:, idx_state(n)));

        xx = ALL.(set_names(n)).GRAND.fly_stats(v).mean.State.mean(:, idx_state(n));
        xfit = xx(fit_range);
        tfit = tt(fit_range);
        [out] = fit_sine(xfit, tfit, freq(v), false);
        %plot(ALL.(set_names(n)).FUNC{v}.All.time, ALL.(set_names(n)).FUNC{v}.All.X, 'Color', [0 0 0 0.5], 'LineWidth', 1)

        FRF_fit.(set_names(n)).(state_names(n)).gain(v,1) = out.mag ./ ALL.(set_names(n)).FUNC{v}.All.Amp;
        FRF_fit.(set_names(n)).(state_names(n)).phase(v,1) = rad2deg(wrapTo2Pi(out.phase - ALL.(set_names(n)).FUNC{v}.All.Phase));
        FRF_fit.(set_names(n)).(state_names(n)).dc(v,1) = out.dc;
        
        if (FRF_fit.(set_names(n)).(state_names(n)).phase(v,1) > 120)
            FRF_fit.(set_names(n)).(state_names(n)).phase(v,1) = FRF_fit.(set_names(n)).(state_names(n)).phase(v,1) - 360;
        end

        % Unclip signal and refit
        if strcmp(set_names(n), "BodyFixed") && strcmp(state_names(n), "head") && any((v == [1:2]), 2) % clipping probaby present here
            [xfit, tfit] = select_point_windows(xfit, tfit, 'y');
            [out] = fit_sine(xfit, tfit, freq(v), false);
            plot(out.x, out.fit, 'Color', [0.5 0.5 0.5 0.5])
            xlim([1 7.7])
            ylim(25*[-1 1])
            pause
            close all

            FRF_fit.(set_names(n)).(state_names(n)).unclip.gain(v,1) = out.mag ./ ALL.(set_names(n)).FUNC{v}.All.Amp;
            FRF_fit.(set_names(n)).(state_names(n)).unclip.phase(v,1) = rad2deg(wrapTo2Pi(out.phase - ALL.(set_names(n)).FUNC{v}.All.Phase));
            FRF_fit.(set_names(n)).(state_names(n)).unclip.dc(v,1) = out.dc;
            
        if (FRF_fit.(set_names(n)).(state_names(n)).unclip.phase(v,1) > 120)
            FRF_fit.(set_names(n)).(state_names(n)).unclip.phase(v,1) = FRF_fit.(set_names(n)).(state_names(n)).unclip.phase(v,1) - 360;
        end
            
        else % just use the regular fit because there was no clipping
            FRF_fit.(set_names(n)).(state_names(n)).unclip.gain(v,1) = FRF_fit.(set_names(n)).(state_names(n)).gain(v,1);
            FRF_fit.(set_names(n)).(state_names(n)).unclip.phase(v,1) = FRF_fit.(set_names(n)).(state_names(n)).phase(v,1);
            FRF_fit.(set_names(n)).(state_names(n)).unclip.dc(v,1) = FRF_fit.(set_names(n)).(state_names(n)).dc(v,1);
        end
        %if FRF_fit.(set_names(n)).(state_names(n)).phase(v,1) > 100
    end
end

%% Comparison
cc.raw = [0.4 0.4 0.6];
cc.lssa = 'r';
cc.unclipped = 'k';

set_names = ["HeadFree", "HeadFree", "BodyFixed"];
state_names = ["body", "head", "head"];
n_set = length(state_names);

close all
fig = figure (406);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_set*2.2 4])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_set);
for n = 1:n_set
    ax(1,n) = subplot(2,n_set,n); cla ; hold on ; title(state_names(n))
        [h.patch(1,n,1),h.line(1,n,1)] = PlotPatch(...
                fft_frf.(set_names(n)).(state_names(n)).gain.mean,...
                fft_frf.(set_names(n)).(state_names(n)).gain.std, ...
                freq, 1, 1, cc.raw, 0.7*cc.raw, 0.2, 1);
        h.line(2,n,1) = plot(freq, FRF_fit.(set_names(n)).(state_names(n)).unclip.gain, 'Color', cc.unclipped);
        h.line(3,n,1) = plot(freq, FRF_fit.(set_names(n)).(state_names(n)).gain, 'Color', cc.lssa);
        if n == 1
            ylabel('gain')
        end

    ax(2,n) = subplot(2,n_set,n+n_set); cla ; hold on
        [h.patch(1,n,2),h.line(1,n,2)] = PlotPatch(...
                fft_frf.(set_names(n)).(state_names(n)).phase.mean,...
                fft_frf.(set_names(n)).(state_names(n)).phase.std, ...
                freq, 1, 1, cc.raw, 0.7*cc.raw, 0.2, 1);
        h.line(2,n,2) = plot(freq, FRF_fit.(set_names(n)).(state_names(n)).unclip.phase, 'Color', cc.unclipped);
        h.line(3,n,2) = plot(freq, FRF_fit.(set_names(n)).(state_names(n)).phase , 'Color', cc.lssa);
        if n == 1
            ylabel('phase difference')
        end
        xlabel('frequency (hz)')
end
% h.line(4) = plot(freq, FRF_fit.BodyFixed.unclip.gain - FRF_fit.BodyFixed.gain, 'k');
linkaxes(ax, 'x')
set(ax, 'xscale', 'log', 'XLim', [0.5 12], 'XTick', [1 10])
set(ax(1,:), 'YLim', [0 1],  'YTick', 0:0.2:1)
set(ax(2,:), 'YLim', [-200 100])
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(h.line(:,:,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)

% xlabel('frequency (hz)')
% ylabel('gain')

leg = legend(h.line(:,1), 'raw', 'unclipped', 'lssa', ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position = [0.599,0.447,0.162,0.138];

set(ax(:,2:end), 'YTickLabels', {})
set(ax(1,:), 'XTickLabels', {})
% set(ax(1,:), 'YScale', 'log')

%% Save FRF data
fname = 'FRF_clipped';
for n = 1:length(FILE)
    filedata = textscan(FILE{n}, '%s', 'delimiter', '._');
    temp_name = [];
    for k = 2:5
        temp_name = [temp_name '_' char(filedata{1}(k))];
    end
    fname = [fname temp_name];
end
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'FRF_fit', 'fft_frf', 'freq');

end