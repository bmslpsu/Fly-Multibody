function [] = SOS_fit_switch()
%% SOS_fit_switch:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

%% Set initialization parameters
clearvars -except PATH FILE ALL MODEL
% sys = MODEL.HeadFree.G.head;
% sys
clss = 'HeadFree';

MODEL.(clss).G.gaze = MODEL.(clss).G.body + MODEL.(clss).G.head;
SWITCH = [];

SWITCH.(clss).G.body = MODEL.HeadFree.P_norm.body * MODEL.HeadFree.C.head;
SWITCH.(clss).G.body.IODelay = MODEL.HeadFree.P_norm.body.IODelay;
SWITCH.(clss).G.head = MODEL.HeadFree.P_norm.head * MODEL.HeadFree.C.body;
SWITCH.(clss).G.head.IODelay = MODEL.HeadFree.G.head.IODelay;
SWITCH.(clss).G.gaze = SWITCH.(clss).G.body + SWITCH.(clss).G.head;

SWITCH.(clss).H.body = minreal( SWITCH.(clss).G.body / (1 + SWITCH.(clss).G.body + SWITCH.(clss).G.head) );
SWITCH.(clss).H.head = minreal( SWITCH.(clss).G.head / (1 + SWITCH.(clss).G.body + SWITCH.(clss).G.head) );
SWITCH.(clss).H.gaze = minreal( (SWITCH.(clss).G.body + SWITCH.(clss).G.head) / (1 + SWITCH.(clss).G.body + SWITCH.(clss).G.head) );

%% Switched closed-loop models comparison
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree", "HeadFree"];
trf = ["ref2body", "ref2head", "ref2gaze"];
cl = ["body", "head", "gaze"];
cc_fit = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1];
cc_scale = [1 0.5];
linespec = ["-", "-"];
n_plot = 1;
n_curve = length(cl);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = 1.5*[n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_curve-1
    sI = 1:4;
  	D = MODEL.(clss(m)).data.(trf(m));
    G = cell(2,1);
    G{1} = (MODEL.(clss(m)).G.(cl(m)));
    G{2} = (SWITCH.(clss(m)).G.(cl(m)));
    %f_cut = G_true.denominator{1}(2) / (2*pi);
    for g = 1:length(G)
        M = plotFit(D.input, IOFv, G{g}, 0:0.02:20, false);
        ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
            title('Nyquist')
            xline(0, '--k')
            yline(0, '--k')
            h.sys(1,m,g) = plot(M(n).real, M(n).imag, 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(1,m,g) = plot(M(n).real_freq, M(n).imag_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)

        ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
            h.sys(2,m,g) = plot(M(n).fv_bode, M(n).gain, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(2,m,g) = plot(D.freq, M(n).gain_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            %xline(f_cut, 'Color', cc_fit(m,:))

        ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
            yline(0, '--k')
            h.sys(3,m,g) = plot(M(n).fv_bode, M(n).phase, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(3,m,g) = plot(D.freq, M(n).phase_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            %xline(f_cut, 'Color', cc_fit(m,:))

        ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
            yline(1, '--k')
            h.sys(4,m,g) = plot(M(n).fv_bode, M(n).error, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(4,m,g) = plot(D.freq, M(n).error_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            ax(4).YLim(1) = 0;
    end
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
% set(ax(2,:), 'YLim', [0 1.1])
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

%% Plant models comparison
clc
clear ax h
IOFv = MODEL.HeadFree.data.err2body.freq;

clss = ["HeadFree", "HeadFree", "HeadFree"];
trf = ["ref2body", "ref2head", "ref2gaze"];
cl = ["body", "head", "gaze"];
cc_fit = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1];
cc_scale = [1 0.5];
linespec = ["-", "-"];
n_plot = 1;
n_curve = length(cl);

fig = figure (2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = 1.5*[n_plot*(3/2) 5];
movegui(fig, 'center')
ax = gobjects(4,n_plot);
n = 1;
for m = 1:n_curve-1
    sI = 1:4;
  	D = MODEL.(clss(m)).data.(trf(m));
    G = cell(1,1);
    G{1} = (MODEL.(clss(m)).P_norm.(cl(m)));
    %f_cut = G_true.denominator{1}(2) / (2*pi);
    for g = 1:length(G)
        M = plotFit(D.input, IOFv, G{g}, 0:0.02:20, false);
        ax(1,m) = subplot(4,n_plot,sI(1)) ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
            title('Nyquist')
            xline(0, '--k')
            yline(0, '--k')
            h.sys(1,m,g) = plot(M(n).real, M(n).imag, 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(1,m,g) = plot(M(n).real_freq, M(n).imag_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)

        ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
            h.sys(2,m,g) = plot(M(n).fv_bode, M(n).gain, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(2,m,g) = plot(D.freq, M(n).gain_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            %xline(f_cut, 'Color', cc_fit(m,:))

        ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
            yline(0, '--k')
            h.sys(3,m,g) = plot(M(n).fv_bode, M(n).phase, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(3,m,g) = plot(D.freq, M(n).phase_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            %xline(f_cut, 'Color', cc_fit(m,:))

        ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
            yline(1, '--k')
            h.sys(4,m,g) = plot(M(n).fv_bode, M(n).error, linespec(g), 'Color', cc_scale(g)*cc_fit(m,:));
            h.sys_freq(4,m,g) = plot(D.freq, M(n).error_freq, ...
                '.', 'Color', cc_scale(g)*cc_fit(m,:), 'MarkerSize', 15);
            ax(4).YLim(1) = 0;
    end
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(2:end,:), 'XLim', [0.1 20])
set(ax(2:end,:), 'XScale', 'log')
set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax(1,:), 'YLim', [-1 1], 'XLim', [-1 1])
% set(ax(2,:), 'YLim', [0 1.1])
set(ax(3,:), 'YLim', [-100 150])
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

end