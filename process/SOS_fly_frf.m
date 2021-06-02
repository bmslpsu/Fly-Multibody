function [] = SOS_fly_frf()
%% SOS_fly_mean:
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
        pI = [1 3 1 2 3 5 6];
        T = ["ref", "ref_wing", "body", "head", "gaze", "dwba", "body_wing"];
        F = ["ref2body","ref2head","ref2gaze","ref2wing","wing2body"];
    case 'BodyFixed'
        pI = [1 1 2 4 5];
        T = ["ref" "head", "dwba", "lwing", "rwing"];
        F = ["ref2head","ref2wing","head2wing"];
end

% Collect fly means
n_state = length(pI);
FLY_mean = [];
FLY_mean.time = DATA.reference{1}.time;
n_detrend = 2;
for n = 1:n_state
    for v = 1:N.vel
        if any(strcmp(T(n), ["ref", "ref_wing"]))
            FLY_mean.fly(v).(T(n)) = squeeze( GRAND.fly_all(v).mean.refState(:,pI(n),:) );
            FLY_mean.grand_mean(v).(T(n)) = GRAND.fly_stats(v).mean.refState.mean(:,pI(n));
            FLY_mean.grand_std(v).(T(n)) = GRAND.fly_stats(v).mean.refState.std(:,pI(n));
        else
            FLY_mean.fly(v).(T(n)) = detrend( squeeze( GRAND.fly_all(v).mean.State(:,pI(n),:)), n_detrend);
            FLY_mean.grand_mean(v).(T(n)) = detrend( GRAND.fly_stats(v).mean.State.mean(:,pI(n)), n_detrend);
            FLY_mean.grand_std(v).(T(n)) = std(FLY_mean.fly(v).(T(n)), [], 2);
        end
    end
end

% Compuete FRF's
FRF_data = [];
n_trf = length(F);
fnames = ["IOMag","IOGain","IOPhaseDiff","IOFRF_error","IOCohr","IOFRF"];
for v = 1:N.vel
    FRF_data.IOFv{v} = GRAND.fly_stats(v).mean.IOFv.mean(1:9);
    for f = 1:size(FLY_mean.fly(v).(T(1)),2)
        % FRF's
        switch clss
            case 'HeadFree'
                FRF_data.fly(v).(F(1))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), ....
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).body(:,f));
                FRF_data.fly(v).(F(2))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).head(:,f));
                FRF_data.fly(v).(F(3))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).gaze(:,f));
                FRF_data.fly(v).(F(4))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref_wing(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).dwba(:,f));
                FRF_data.fly(v).(F(5))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).dwba(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).body_wing(:,f));
            case 'BodyFixed'
                FRF_data.fly(v).(F(1))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).head(:,f));
                FRF_data.fly(v).(F(2))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).ref(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).dwba(:,f));
                FRF_data.fly(v).(F(3))(f,1) = frf(FLY_mean.time, FLY_mean.fly(v).head(:,f), ...
                    FRF_data.IOFv{v}, false, FLY_mean.fly(v).dwba(:,f));  
        end

        % Store mag, gain, phase, coherence
        for k = 1:n_trf
            for j = 1:length(fnames)
                if strcmp(fnames(j),"IOPhaseDiff")
                    phs = rad2deg(FRF_data.fly(v).(F(k))(f,1).(fnames(j)));
                    if strcmp(F(k),"ref2body")
                        phs(phs > 30) = phs(phs > 30) - 360;
                    end
                    if strcmp(F(k),"ref2head")
                        %phs(phs < -100) = phs(phs < -100) + 360;
                        phs(phs < -170) = phs(phs < -170) + 360;
                    end
                    if strcmp(F(k),"wing2body")
                        phs(phs < -50) = phs(phs < -50) + 360;
                        %phs(phs > 0) = phs(phs > 0) - 360;
                        %phs = phs + 0;
                    end
                    if strcmp(F(k),"ref2wing")
                        %phs(phs > 10) = phs(phs > 10) - 360;
                    end
                    FRF_data.fly_all.(F(k))(v).(fnames(j))(:,f) = phs;
                else
                    FRF_data.fly_all.(F(k))(v).(fnames(j))(:,f) = FRF_data.fly(v).(F(k))(f,1).(fnames(j));
                end
                
                cohr_rmv =  FRF_data.fly(v).(F(k))(f,1).IOCohr < 0.6;
                FRF_data.fly_all.(F(k))(v).(fnames(j))(cohr_rmv,f) = nan;
                    
                FRF_data.grand_mean.(F(k))(v).(fnames(j)) = nanmean(FRF_data.fly_all.(F(k))(v).(fnames(j)), 2);
                FRF_data.grand_std.(F(k))(v).(fnames(j)) = nanstd(FRF_data.fly_all.(F(k))(v).(fnames(j)), [], 2);
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
v = 1;
for n = 1:n_trf
    subI = n + (0:4)*n_trf;
    ax(1,n) = subplot(5,n_trf,subI(1)); hold on ; title(F(n), 'interpreter', 'none')
        plot(FRF_data.IOFv{v}, FRF_data.fly_all.(F(n))(v).IOMag, '.', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(1,n),h.line(1,n)] = PlotPatch(FRF_data.grand_mean.(F(n))(v).IOMag,...
                  FRF_data.grand_std.(F(n))(v).IOMag, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(2,n) = subplot(5,n_trf,subI(2)); hold on
        plot(FRF_data.IOFv{v}, FRF_data.fly_all.(F(n))(v).IOGain, '.', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(2,n),h.line(2,n)] = PlotPatch(FRF_data.grand_mean.(F(n))(v).IOGain,...
                  FRF_data.grand_std.(F(n))(v).IOGain, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(3,n) = subplot(5,n_trf,subI(3)); hold on
        plot([0 20], [0 0], '--k')
        plot(FRF_data.IOFv{v}, FRF_data.fly_all.(F(n))(v).IOPhaseDiff, '.', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(3,n),h.line(3,n)] = PlotPatch(FRF_data.grand_mean.(F(n))(v).IOPhaseDiff,...
                  FRF_data.grand_std.(F(n))(v).IOPhaseDiff, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(4,n) = subplot(5,n_trf,subI(4)); hold on
        plot(FRF_data.IOFv{v}, FRF_data.fly_all.(F(n))(v).IOFRF_error, '.', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.grand_mean.(F(n))(v).IOFRF_error,...
                  FRF_data.grand_std.(F(n))(v).IOFRF_error, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
    ax(5,n) = subplot(5,n_trf,subI(5)); hold on
        plot(FRF_data.IOFv{v}, FRF_data.fly_all.(F(n))(v).IOCohr, '.', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.grand_mean.(F(n))(v).IOCohr,...
                  FRF_data.grand_std.(F(n))(v).IOCohr, FRF_data.IOFv{v}, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
end

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15')
set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.2 25],...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
% linkaxes(ax(1,:),'y')
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
set(ax(3,1:end),'YLim',[-300 300])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])

set(ax,'XScale','log')
align_Ylabels(fig)

%% Save Fly mean & FRF data
filedata = textscan(FILE, '%s', 'delimiter', '_');
dataset_name = [];
for n = 1:6
    dataset_name = [dataset_name '_' char(filedata{1}(n))];
end
fname = ['FRF_fly_mean' dataset_name];
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'FLY_mean', 'FRF_data', 'FUNC', 'U', 'N');
end