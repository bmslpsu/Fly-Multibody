function [] = SS_wing2body_bode()
%% SS_wing2body_bode:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat'}, ...
    'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','FUNC','GRAND','FLY','D','I','U','N')

%% Collect fly means
close all
clc
clearvars -except FILE DATA ALL GRAND FLY FUNC D I U N root
% [Mag,Phase] = computeMagPhase(time, meanData, f)
time = GRAND.fly_all(1).mean.Time(:,1);
Fs = 1 / mean(diff(time));
Fv = GRAND.fly_all(1).mean.Fv(:,1);
IOFv = U.freq{1};

pI = [1 2 3 5];
T = ["body", "head", "gaze", "wing", "ref"];
ALL = [];
for n = 1:length(T)
    for v = 1:N.freq
        if n == length(T)
            ALL.fly(v).(T(n)) = squeeze( GRAND.fly_all(v).mean.refState(:,1,:) );
            ALL.grand_mean(v).(T(n)) = GRAND.fly_stats(v).mean.refState.mean(:,1);
            ALL.grand_std(v).(T(n)) = GRAND.fly_stats(v).mean.refState.std(:,1);
        else
            ALL.fly(v).(T(n)) = squeeze( GRAND.fly_all(v).mean.State(:,pI(n),:) );
            ALL.grand_mean(v).(T(n)) = GRAND.fly_stats(v).mean.State.mean(:,pI(n));
            ALL.grand_std(v).(T(n)) = GRAND.fly_stats(v).mean.State.std(:,pI(n));
        end        
        %freq_field = [char(T(n)) '_freq'];
        %mag_field = [char(T(n)) '_mag'];
        %ALL.grand_mean(v).(freq_field) = chirpz(ALL.grand_mean(v).(T(n)),Fs,0,Fs/2);
    end
end
%% ref2wing
FRF_data = [];
for v = 1:N.freq
    for f = 1:N.fly
        FRF_data.fly(v).wing2body(f) = frf(time, ALL.fly(v).ref(:,f), IOFv(v), false, ...
                                                ALL.fly(v).wing(:,f));
        if FRF_data.fly(v).wing2body(f).IOCohr > 0.5
            FRF_data.fly(v).wing2body_mag(f) = FRF_data.fly(v).wing2body(f).refIOMag;
            FRF_data.fly(v).wing2body_gain(f) = FRF_data.fly(v).wing2body(f).IOGain;
            FRF_data.fly(v).wing2body_cohr(f) = FRF_data.fly(v).wing2body(f).IOCohr;
            FRF_data.fly(v).wing2body_phase(f) = rad2deg(FRF_data.fly(v).wing2body(f).IOPhaseDiff);
        else
            FRF_data.fly(v).wing2body_mag(f) = nan;
            FRF_data.fly(v).wing2body_gain(f) = nan;
            FRF_data.fly(v).wing2body_cohr(f) = nan;
            FRF_data.fly(v).wing2body_phase(f) = nan;
        end
        
%         if FRF_data.fly(v).wing2body_phase(f) < -10
%             FRF_data.fly(v).wing2body_phase(f) = FRF_data.fly(v).wing2body_phase(f) + 360;
%         end
    end
    
    fnames = fieldnames(FRF_data.fly);
    for n = 2:length(fnames)
        FRF_data.fly_all.(fnames{n}) = cat(1, FRF_data.fly(:).(fnames{n}));
        FRF_data.fly_mean.(fnames{n}) = nanmean(FRF_data.fly_all.(fnames{n}), 2);
        FRF_data.fly_std.(fnames{n}) = nanstd(FRF_data.fly_all.(fnames{n}), [], 2);
    end

    FRF_data.grand_mean(v).wing2body = frf(time, ALL.grand_mean(v).ref, IOFv(v), false, ALL.grand_mean(v).wing);
    FRF_data.wing2body_mag(v) = FRF_data.grand_mean(v).wing2body.refIOMag;
    FRF_data.wing2body_gain(v) = FRF_data.grand_mean(v).wing2body.IOGain;
    FRF_data.wing2body_cohr(v) = FRF_data.grand_mean(v).wing2body.IOCohr;
    FRF_data.wing2body_phase(v) = rad2deg(FRF_data.grand_mean(v).wing2body.IOPhaseDiff);
    
%     if FRF_data.wing2body_phase(v) < -20
%         FRF_data.wing2body_phase(v) = FRF_data.wing2body_phase(v) + 360;
%     end
end

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 4*1.8])
movegui(fig, 'center')
clear ax h
cc = [0.9 0.1 0.8];
ax = gobjects(4,1);
ax(1,1) = subplot(4,1,1); hold on
    plot(IOFv, FRF_data.wing2body_mag, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_mag, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(1,1),h.line(1,1)] = PlotPatch(FRF_data.fly_mean.wing2body_mag,...
              FRF_data.fly_std.wing2body_mag, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(2,1) = subplot(4,1,2); hold on
    plot(IOFv, FRF_data.wing2body_gain, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_gain, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(2,1),h.line(2,1)] = PlotPatch(FRF_data.fly_mean.wing2body_gain,...
              FRF_data.fly_std.wing2body_gain, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(3,1) = subplot(4,1,3); hold on
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
    plot(IOFv, FRF_data.wing2body_phase, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_phase, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(3,1),h.line(3,1)] = PlotPatch(FRF_data.fly_mean.wing2body_phase,...
              FRF_data.fly_std.wing2body_phase, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(4,1) = subplot(4,1,4); hold on
    plot(IOFv, FRF_data.wing2body_cohr, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_cohr, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(4,1),h.line(4,1)] = PlotPatch(FRF_data.fly_mean.wing2body_cohr,...
              FRF_data.fly_std.wing2body_cohr, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.5 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
% set(ax, 'XTick', [1:10])

% delete(h.patch)

% set(ax(1), 'YLim', [0 2])
% set(ax(2), 'YLim', [0 100])
% set(ax(3), 'YLim', [0 100])
set(ax(4), 'YLim', [0 1])

set(ax,'XScale','log')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°/s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

set(ax(1:end-1,:), 'XColor', 'none')



%% wing2body
FRF_data = [];
for v = 1:N.freq
    for f = 1:N.fly
        FRF_data.fly(v).wing2body(f) = frf(time, ALL.fly(v).wing(:,f), IOFv(v), false, ...
                                                ALL.fly(v).body(:,f));
        if FRF_data.fly(v).wing2body(f).IOCohr > 0.5
            FRF_data.fly(v).wing2body_mag(f) = FRF_data.fly(v).wing2body(f).refIOMag;
            FRF_data.fly(v).wing2body_gain(f) = FRF_data.fly(v).wing2body(f).IOGain;
            FRF_data.fly(v).wing2body_cohr(f) = FRF_data.fly(v).wing2body(f).IOCohr;
            FRF_data.fly(v).wing2body_phase(f) = rad2deg(FRF_data.fly(v).wing2body(f).IOPhaseDiff);
        else
            FRF_data.fly(v).wing2body_mag(f) = nan;
            FRF_data.fly(v).wing2body_gain(f) = nan;
            FRF_data.fly(v).wing2body_cohr(f) = nan;
            FRF_data.fly(v).wing2body_phase(f) = nan;
        end
        
        if FRF_data.fly(v).wing2body_phase(f) > 150
            FRF_data.fly(v).wing2body_phase(f) = FRF_data.fly(v).wing2body_phase(f) - 360;
        end
    end
    
    fnames = fieldnames(FRF_data.fly);
    for n = 2:length(fnames)
        FRF_data.fly_all.(fnames{n}) = cat(1, FRF_data.fly(:).(fnames{n}));
        FRF_data.fly_mean.(fnames{n}) = nanmean(FRF_data.fly_all.(fnames{n}), 2);
        FRF_data.fly_std.(fnames{n}) = nanstd(FRF_data.fly_all.(fnames{n}), [], 2);
    end

    FRF_data.grand_mean(v).wing2body = frf(time, ALL.grand_mean(v).wing, IOFv(v), false, ALL.grand_mean(v).body);
    FRF_data.wing2body_mag(v) = FRF_data.grand_mean(v).wing2body.refIOMag;
    FRF_data.wing2body_gain(v) = FRF_data.grand_mean(v).wing2body.IOGain;
    FRF_data.wing2body_cohr(v) = FRF_data.grand_mean(v).wing2body.IOCohr;
    FRF_data.wing2body_phase(v) = rad2deg(FRF_data.grand_mean(v).wing2body.IOPhaseDiff);
    
    if FRF_data.wing2body_phase(v) > 150
        FRF_data.wing2body_phase(v) = FRF_data.wing2body_phase(v) - 360;
    end
end

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 4*1.8])
movegui(fig, 'center')
clear ax h
cc = [0.9 0.1 0.8];
ax = gobjects(4,1);
ax(1,1) = subplot(4,1,1); hold on
    plot(IOFv, FRF_data.wing2body_mag, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_mag, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(1,1),h.line(1,1)] = PlotPatch(FRF_data.fly_mean.wing2body_mag,...
              FRF_data.fly_std.wing2body_mag, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(2,1) = subplot(4,1,2); hold on
    plot(IOFv, FRF_data.wing2body_gain, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_gain, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(2,1),h.line(2,1)] = PlotPatch(FRF_data.fly_mean.wing2body_gain,...
              FRF_data.fly_std.wing2body_gain, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(3,1) = subplot(4,1,3); hold on
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
    plot(IOFv, FRF_data.wing2body_phase, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_phase, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(3,1),h.line(3,1)] = PlotPatch(FRF_data.fly_mean.wing2body_phase,...
              FRF_data.fly_std.wing2body_phase, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);
ax(4,1) = subplot(4,1,4); hold on
    plot(IOFv, FRF_data.wing2body_cohr, '-k.', 'LineWidth', 1, 'MarkerSize', 10)
    plot(IOFv, FRF_data.fly_all.wing2body_cohr, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
    [h.patch(4,1),h.line(4,1)] = PlotPatch(FRF_data.fly_mean.wing2body_cohr,...
              FRF_data.fly_std.wing2body_cohr, IOFv, 1, 1, cc, 0.7*cc, 0.2, 1);

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.5 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
% set(ax, 'XTick', [1:10])

% delete(h.patch)

set(ax(1), 'YLim', [0 2])
% set(ax(2), 'YLim', [0 100])
% set(ax(3), 'YLim', [0 100])
set(ax(4), 'YLim', [0 1])

set(ax,'XScale','log')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°/s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

set(ax(1:end-1,:), 'XColor', 'none')


%% Fit data
close all
clc

start_time = 2;
startI = start_time*Fs;
n_detrend = 7;
for n = 1:length(T)
    disp(T(n))
    mag_field = [char(T(n)) '_mag'];
    phase_field = [char(T(n)) '_phase'];
    gof_field = [char(T(n)) '_gof'];
    for v = 1:N.freq
        for f = 1:N.fly
            y = ALL.fly(v).(T(n))(startI:end,f);
            [fitresult, gof] = fitSS(time(startI:end), y, FUNC{v}.All.Freq, n_detrend, false);
            
            if gof.rsquare > 0.5
                ALL.fly(v).(gof_field)(f) = gof.rsquare;
                ALL.fly(v).(mag_field)(f) = fitresult.a;
                ALL.fly(v).(phase_field)(f) = fitresult.b;
                
            else
                ALL.fly(v).(gof_field)(f) = nan;
                ALL.fly(v).(mag_field)(f) = nan;
                ALL.fly(v).(phase_field)(f) = nan;
            end
            %pause
            %cla
        end
    end
end

%% Save FRF data
% filedata = textscan(FILE, '%s', 'delimiter', '_');
% dataset_name = [];
% for n = 1:4
%     dataset_name = [dataset_name '_' char(filedata{1}(n))];
% end
% fname = ['FRF' dataset_name];
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'FRF_data', 'FUNC', 'U', 'N');

end