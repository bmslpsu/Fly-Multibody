function [] = SS_fly_frf()
%% SS_fly_mean:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat'},'Select file', root, 'MultiSelect','off');
filedata = textscan(FILE, '%s', 'delimiter', '_');
clss = filedata{1}{2};
load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N')

%% System ID on fly means
clearvars -except FUNC DATA ALL GRAND FLY D I U N root FILE clss
clc

switch clss
    case 'HeadFree'
        pI = [1 1 2 3 5 10 11];
        T = ["ref", "body", "head", "gaze", "dwba", "lwing", "rwing"];
        F = ["ref2body","ref2head","ref2gaze","ref2wing","wing2body"]; % FRF's
    case 'BodyFixed'
        pI = [1 1 2 4 5];
        T = ["ref" "head", "dwba", "lwing", "rwing"];
        F = ["ref2head","ref2wing","head2wing"]; % FRF's
end

% Collect fly means
n_state = length(pI);
FLY_mean = [];
FLY_mean.time = DATA.reference{1}.time;
n_detrend = 2;
for n = 1:n_state
    for v = 1:N.freq
        if strcmp(T(n), 'ref')
            FLY_mean.fly(v).(T(n)) = squeeze( GRAND.fly_all(v).mean.refState(:,1,:) );
            FLY_mean.grand_mean(v).(T(n)) = GRAND.fly_stats(v).mean.refState.mean(:,1);
            FLY_mean.grand_std(v).(T(n)) = GRAND.fly_stats(v).mean.refState.std(:,1);
        else
            FLY_mean.fly(v).(T(n)) = detrend( squeeze( GRAND.fly_all(v).mean.State(:,pI(n),:)), n_detrend);
            FLY_mean.grand_mean(v).(T(n)) = detrend( GRAND.fly_stats(v).mean.State.mean(:,pI(n)), n_detrend);
            FLY_mean.grand_std(v).(T(n)) = std(FLY_mean.fly(v).(T(n)), [], 2);
        end
    end
end

% Compuete FRF's
FRF_data = [];
FRF_data.IOFv = U.freq{1};
n_trf = length(F);
fnames = ["IOMag","IOGain","IOPhaseDiff","IOFRF_error","IOCohr","IOFRF"];
for v = 1:N.freq
    for f = 1:size(FLY_mean.fly(v).(T(1)),2)
        freq = FUNC{v}.All.Freq;

        % FRF's
        switch clss
            case 'HeadFree'
                FRF_data.fly(v).(F(1))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).body(:,f));
                FRF_data.fly(v).(F(2))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).head(:,f));
                FRF_data.fly(v).(F(3))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).gaze(:,f));
                FRF_data.fly(v).(F(4))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).dwba(:,f));
                FRF_data.fly(v).(F(5))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).dwba(:,f), freq, false, FLY_mean.fly(v).body(:,f));
            case 'BodyFixed'
                FRF_data.fly(v).(F(1))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).head(:,f));
                FRF_data.fly(v).(F(2))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), freq, false, FLY_mean.fly(v).dwba(:,f));
                FRF_data.fly(v).(F(3))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).head(:,f), freq, false, FLY_mean.fly(v).dwba(:,f));  
        end

        % Store mag, gain, phase, coherence
        for k = 1:n_trf
            for j = 1:length(fnames)
                if  FRF_data.fly(v).(F(k))(f,1).IOCohr > 0.6
                    if strcmp(fnames(j),"IOPhaseDiff")
                        phs = rad2deg(FRF_data.fly(v).(F(k))(f,1).(fnames(j)));
                        if strcmp(F(k),"wing2body")
                            phs(phs < -50) = phs(phs < -50) + 360;
                        end
                        if strcmp(F(k),"ref2wing")
                            phs(phs > 0) = phs(phs > 0) - 360;
                        end
                        FRF_data.fly_all.(F(k)).(fnames(j))(v,f) = phs;
                    else
                        FRF_data.fly_all.(F(k)).(fnames(j))(v,f) = FRF_data.fly(v).(F(k))(f,1).(fnames(j));
                    end
                else
                    FRF_data.fly_all.(F(k)).(fnames(j))(v,f) = nan;
                end
                FRF_data.grand_mean.(F(k)).(fnames(j)) = nanmean(FRF_data.fly_all.(F(k)).(fnames(j)), 2);
                FRF_data.grand_std.(F(k)).(fnames(j)) = nanstd(FRF_data.fly_all.(F(k)).(fnames(j)), [], 2);
            end
        end
    end
end

%% FRF: one condition
cc = hsv(n_trf);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.2*n_trf 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,n_trf);
for n = 1:n_trf
    subI = n + (0:4)*n_trf;
    ax(1,n) = subplot(5,n_trf,subI(1)); hold on ; title(F(n), 'interpreter', 'none')
        plot(FRF_data.IOFv, FRF_data.fly_all.(F(n)).IOMag, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(1,n),h.line(1,n)] = PlotPatch(FRF_data.grand_mean.(F(n)).IOMag,...
                  FRF_data.grand_std.(F(n)).IOMag, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(2,n) = subplot(5,n_trf,subI(2)); hold on
        plot(FRF_data.IOFv, FRF_data.fly_all.(F(n)).IOGain, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(2,n),h.line(2,n)] = PlotPatch(FRF_data.grand_mean.(F(n)).IOGain,...
                  FRF_data.grand_std.(F(n)).IOGain, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(3,n) = subplot(5,n_trf,subI(3)); hold on
        plot([0 20], [0 0], '--k')
        plot(FRF_data.IOFv, FRF_data.fly_all.(F(n)).IOPhaseDiff, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(3,n),h.line(3,n)] = PlotPatch(FRF_data.grand_mean.(F(n)).IOPhaseDiff,...
                  FRF_data.grand_std.(F(n)).IOPhaseDiff, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(4,n) = subplot(5,n_trf,subI(4)); hold on
        plot(FRF_data.IOFv, FRF_data.fly_all.(F(n)).IOFRF_error, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.grand_mean.(F(n)).IOFRF_error,...
                  FRF_data.grand_std.(F(n)).IOFRF_error, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
    ax(5,n) = subplot(5,n_trf,subI(5)); hold on
        plot(FRF_data.IOFv, FRF_data.fly_all.(F(n)).IOCohr, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.grand_mean.(F(n)).IOCohr,...
                  FRF_data.grand_std.(F(n)).IOCohr, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
end

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15')
set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.5 25],...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(3,:), 'y')
linkaxes(ax(4,:), 'y')
linkaxes(ax(5,:), 'y')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

set(ax(2,1:end-1),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-250 200])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])

set(ax,'XScale','log')
align_Ylabels(fig)

%% Save Fly mean & FRF data
filedata = textscan(FILE, '%s', 'delimiter', '_');
dataset_name = [];
for n = 1:5
    dataset_name = [dataset_name '_' char(filedata{1}(n))];
end
fname = ['FRF_fly_mean' dataset_name];
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'FLY_mean', 'FRF_data', 'FUNC', 'U', 'N');
end