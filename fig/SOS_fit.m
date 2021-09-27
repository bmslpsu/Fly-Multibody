function [] = SOS_fit()
%% SOS_fit:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

%% Set initialization parameters
clearvars -except PATH FILE ALL MODEL
% sys = MODEL.HeadFree.G.head;
% sys

%% err2body & err2head with data head-free
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;
n_freq = length(IOFv);

clss = ["HeadFree", "HeadFree"];
trf = ["err2body", "err2head"];
% trf = ["ref2body", "ref2head"];
cc_fit = [0.9 0 0 ; 0 0.4 1];
cc_data = repmat([0.5 0.5 0.5], [n_freq, 1]);
n_plot = length(trf);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:3)*n_plot;
    M = MODEL.(clss(m)).fit.(trf(m));
    D = MODEL.(clss(m)).data.(trf(m));
    ax(1,m) = subplot(4,n_plot,sI(1)); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        h.data(1,m,:) = gscatter(D.real(:), D.imag(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)); cla ; hold on ; box on
        h.data(2,m,:) = gscatter(D.freq_all(:), D.gain(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)); cla ; hold on ; box on
        yline(0, '--k')
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        yline(1, '--k')
        h.data(4,m,:) = gscatter(D.freq_all(:), D.error(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(2,1), 'YLim', [0 3])
set(ax(2,2), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% err2body & err2head with data head-fixed & body_fixed
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;
n_freq = length(IOFv);

clss = ["HeadFixed", "BodyFixed"];
trf = ["err2body", "err2head"];
cc_fit = [1 0.6 0.1 ; 0 0.8 0.2];
cc_data = repmat([0.5 0.5 0.5], [n_freq, 1]);
n_plot = length(trf);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:3)*n_plot;
    M = MODEL.(clss(m)).fit.(trf(m));
    D = MODEL.(clss(m)).data.(trf(m));
    ax(1,m) = subplot(4,n_plot,sI(1)); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        h.data(1,m,:) = gscatter(D.real(:), D.imag(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)); cla ; hold on ; box on
        h.data(2,m,:) = gscatter(D.freq_all(:), D.gain(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)); cla ; hold on ; box on
        yline(0, '--k')
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        yline(1, '--k')
        h.data(4,m,:) = gscatter(D.freq_all(:), D.error(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(2,1), 'YLim', [0 3])
set(ax(2,2), 'YLim', [0 1.2])
set(ax(3,:), 'YLim', [-250 150])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% ref2body & ref2head with data from err2body & err2head
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;
n_freq = length(IOFv);

clss = ["HeadFree", "HeadFree", "HeadFree"];
trf = ["ref2body", "ref2head", "ref2gaze"];
cl = ["body", "head", "gaze"];
cc_fit = [0.9 0 0 ; 0 0.4 1;  0.5 0.3 1];
cc_data = repmat([0.5 0.5 0.5], [n_freq, 1]);
n_plot = length(trf);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:3)*n_plot;
    D = MODEL.(clss(m)).data.(trf(m));
	M = plotFit(D.input, IOFv, MODEL.(clss(m)).H.(cl(m)), 0:0.02:20, false);
    ax(1,m) = subplot(4,n_plot,sI(1)); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k');
        yline(0, '--k');
        h.data(1,m,:) = gscatter(D.real(:), D.imag(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)); cla ; hold on ; box on
        h.data(2,m,:) = gscatter(D.freq_all(:), D.gain(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)); cla ; hold on ; box on
        yline(0, '--k');
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        yline(1, '--k');
        h.data(4,m,:) = gscatter(D.freq_all(:), D.error(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
set(ax(2,:), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabel', [])
set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% ref2body & ref2head with data from err2body & err2head head-fixed & body-fixed
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;
n_freq = length(IOFv);

clss = ["HeadFixed", "BodyFixed"];
trf = ["ref2body", "ref2head"];
cl = ["body", "head"];
cc_fit = [1 0.6 0.1 ; 0 0.8 0.2];
cc_data = repmat([0.5 0.5 0.5], [n_freq, 1]);
n_plot = length(trf);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:3)*n_plot;
    D = MODEL.(clss(m)).data.(trf(m));
	M = plotFit(D.input, IOFv, MODEL.(clss(m)).H.(cl(m)), 0:0.02:20, false);
    ax(1,m) = subplot(4,n_plot,sI(1)); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        h.data(1,m,:) = gscatter(D.real(:), D.imag(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)); cla ; hold on ; box on
        h.data(2,m,:) = gscatter(D.freq_all(:), D.gain(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)); cla ; hold on ; box on
        yline(0, '--k')
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        yline(1, '--k')
        h.data(4,m,:) = gscatter(D.freq_all(:), D.error(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
set(ax(2,:), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabel', [])
set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% ref2body for head-free and head-fixed
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree", "HeadFixed"];
trf = ["ref2body", "ref2head", "ref2body"];
cl = ["body", "body_no_head", "body"];
cc_fit = [0.9 0 0 ; 0 0 0 ;  1 0.6 0.1];
n_plot = 1;
n_curve = length(cl);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_curve
    sI = 1:4;
    D = MODEL.(clss(m)).data.(trf(m));
	M = plotFit(D.input, IOFv, MODEL.(clss(m)).H.(cl(m)), 0:0.02:20, false);
    ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
        yline(0, '--k')
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
        yline(1, '--k')
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
set(ax(2,:), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabels', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% ref2head for head-free and head-fixed
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree", "BodyFixed"];
trf = ["ref2head", "ref2head", "ref2head"];
cl = ["head", "head_no_body", "head"];
cc_fit = [0 0.4 1 ; 0 0 0 ; 0 0.8 0.2];
n_plot = 1;
n_curve = length(cl);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_curve
    sI = 1:4;
    D = MODEL.(clss(m)).data.(trf(m));
	M = plotFit(D.input, IOFv, MODEL.(clss(m)).H.(cl(m)), 0:0.02:20, false);
    ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
        yline(0, '--k')
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
        yline(1, '--k')
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
set(ax(2,:), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabels', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% Head & body controller comparison
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree"];
trf = ["ref2body", "ref2head"];
cl = ["body", "head"];
cc_fit = [0.9 0 0 ; 0 0.4 1];
n_plot = 1;
n_curve = length(cl);

fig = figure (30); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_curve
    sI = 1:4;
  	D = MODEL.(clss(m)).data.(trf(m));
	M = plotFit(D.input, IOFv, MODEL.(clss(m)).C_norm.(cl(m)), 0:0.02:20, false);
    ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k');
        yline(0, '--k');
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
        yline(0, '--k');
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
        yline(1, '--k');
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
% set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
% set(ax(2,:), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-100 100])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabels', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% Head & body plant comparison
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree"];
trf = ["ref2body", "ref2head"];
cl = ["body", "head"];
cc_fit = [0.9 0 0 ; 0 0.4 1];
n_plot = 1;
n_curve = length(cl);

fig = figure (2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
temp = [];
for m = 1:n_curve
    sI = 1:4;
  	D = MODEL.(clss(m)).data.(trf(m));
    P = MODEL.(clss(m)).P_norm.(cl(m));
    f_cut = (1 / P.denominator{1}(1)) / (2*pi);
	M = plotFit(D.input, IOFv, P, 0:0.02:20, false);
    %temp(:,m) = M(n).gain;
    ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k');
        yline(0, '--k');
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        xline(f_cut, 'Color', cc_fit(m,:));
        
    ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
        yline(0, '--k');
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        xline(f_cut, 'Color', cc_fit(m,:));

    ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
        yline(1, '--k');
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
% set(ax(2,:), 'YLim', [0 1.2])
set(ax(3,:), 'YLim', [-100 100])
set(ax(4,:), 'YLim', [0 1.5])
linkaxes(ax(2:end,:), 'x')
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabels', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

%% err2head with proportional controller
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;
n_freq = length(IOFv);

clss = ["HeadFree", "HeadFree"];
trf = ["err2body", "err2head"];
cc_fit = [0.9 0 0 ; 0 0.4 1];
cc_data = repmat([0.5 0.5 0.5], [n_freq, 1]);
n_plot = length(trf);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 2;
for m = 1:n_plot
    sI = m + (0:3)*n_plot;
    M = MODEL.(clss(m)).fit.(trf(m));
    D = MODEL.(clss(m)).data.(trf(m));
    ax(1,m) = subplot(4,n_plot,sI(1)); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k');
        yline(0, '--k');
        h.data(1,m,:) = gscatter(D.real(:), D.imag(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)); cla ; hold on ; box on
        h.data(2,m,:) = gscatter(D.freq_all(:), D.gain(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)); cla ; hold on ; box on
        yline(0, '--k');
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        yline(1, '--k');
        h.data(4,m,:) = gscatter(D.freq_all(:), D.error(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(4,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit(m,:));
        h.sys_freq(4,m) = plot(D.freq, M(n).error_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        ax(4).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(2,1), 'YLim', [0 3])
set(ax(2,2), 'YLim', [0 1])
set(ax(3,:), 'YLim', [-250 150])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6)

set(ax(2:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Real')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Imaginary')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')
XLabelHC = get(ax(4,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

end