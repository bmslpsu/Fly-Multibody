function [] = SOS_fit_v2()
%% SOS_fit_v2:
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

clss = ["HeadFree", "HeadFree", "BodyFixed"];
trf = ["err2body", "err2head", "err2head"];
cc_fit = [0.9 0 0 ; 0 0.4 1 ; 0 0.8 0.2];
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
        %xline(0, '--k');
        %yline(0, '--k');
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
        %yline(0, '--k');
        h.data(3,m,:) = gscatter(D.freq_all(:), D.phase(:), D.freq_all(:), ...
            cc_data, '.', 10, false);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)); cla ; hold on ; box on
        %yline(1, '--k');
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
set(ax(2,2:3), 'YLim', [0 1.2])
set(ax(3,:), 'YLim', [-250 150])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 4)
set(h.sys_freq, 'MarkerSize', 6, 'LineWidth', 0.75)
set(h.sys_freq, 'Marker', 'none')

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

fig = figure (2); clf
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
        %xline(0, '--k');
        %yline(0, '--k');
        h.sys(1,m) = plot(M(n).real, M(n).imag, 'Color', cc_fit(m,:));
        h.sys_freq(1,m) = plot(M(n).real_freq, M(n).imag_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        plot(M(n).real(1), M(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        
    ax(2,m) = subplot(4,n_plot,sI(2)) ; hold on ; box on
        h.sys(2,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit(m,:));
        h.sys_freq(2,m) = plot(D.freq, M(n).gain_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);
        
    ax(3,m) = subplot(4,n_plot,sI(3)) ; hold on ; box on
        %yline(0, '--k');
        h.sys(3,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit(m,:));
        h.sys_freq(3,m) = plot(D.freq, M(n).phase_freq, ...
            '.', 'Color', cc_fit(m,:), 'MarkerSize', 15);

    ax(4,m) = subplot(4,n_plot,sI(4)) ; hold on ; box on
        %yline(1, '--k')
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
set(h.sys_freq, 'MarkerSize', 6, 'LineWidth', 0.75)
set(h.sys_freq, 'Marker', 'none')

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

%% Adaptive head control
G_head_fixed_no_delay = MODEL.BodyFixed.G.head;
G_head_no_delay = MODEL.HeadFree.G.head;

delay_diff = G_head_fixed_no_delay.IODelay - G_head_no_delay.IODelay;

G_head_fixed_no_delay.IODelay = 0;
G_head_no_delay.IODelay = 0;

MODEL.BodyFixed.C_adaptive = G_head_fixed_no_delay / G_head_no_delay;
MODEL.BodyFixed.C_adaptive.IODelay = delay_diff;

for v = 1:3
    ALL.BodyFixed.FRF_data.C_adaptive(v).complex = ALL.BodyFixed.FRF_data.err2head.grand_mean(v).complex ./ ...
                                                    ALL.HeadFree.FRF_data.err2head.grand_mean(v).complex;
                                                
	ALL.BodyFixed.FRF_data.C_adaptive(v).gain = abs(ALL.BodyFixed.FRF_data.C_adaptive(v).complex);
    ALL.BodyFixed.FRF_data.C_adaptive(v).phase = rad2deg( angle(ALL.BodyFixed.FRF_data.C_adaptive(v).complex) );
    ALL.BodyFixed.FRF_data.C_adaptive(v).error = abs( (1 + 0*1i) - ALL.BodyFixed.FRF_data.C_adaptive(v).complex);
end

vI = 2;
% gain_ratio = ALL.BodyFixed.FRF_data.err2head.grand_mean(vI).gain ...
%             ./ ALL.HeadFree.FRF_data.err2head.grand_mean(vI).gain;
% phase_diff = ALL.BodyFixed.FRF_data.err2head.grand_mean(vI).phase ...
%             - ALL.HeadFree.FRF_data.err2head.grand_mean(vI).phase;
% D = {gain_ratio, phase_diff};
D = {ALL.BodyFixed.FRF_data.C_adaptive(vI).gain, ALL.BodyFixed.FRF_data.C_adaptive(vI).phase};
M = plotFit(D, IOFv, MODEL.BodyFixed.C_adaptive, 0:0.02:20, false);

% data = frd(ALL.BodyFixed.FRF_data.C_adaptive(vI).complex, IOFv, 'FrequencyUnit', 'hz');
% opt = tfestOptions('EnforceStability', false, 'InitializeMethod', 'all');
% tffit = [];
% tffit{end+1} = tfest(data, 1, 1, 0.009, opt);
% tffit{1}
% M = plotFit(D, IOFv, tffit{1}, 0:0.02:20, false);

clc
clear ax h

cc_fit = 'c';
cc_data = 'k';
n_plot = 1;

fig = figure (3); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) (5/4)*2.5];
movegui(fig, 'center')
ax = gobjects(3,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:2)*n_plot;        
    ax(1,m) = subplot(3,n_plot,sI(1)); cla ; hold on ; box on
        h.data(1,m,:) = plot(IOFv, D{1}, '.-', 'Color', cc_data, 'MarkerSize', 15);
        h.sys(1,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit);
        
    ax(2,m) = subplot(3,n_plot,sI(2)); cla ; hold on ; box on
        yline(0, '--k');
        h.data(2,m,:) = plot(IOFv, D{2}, '.-', 'Color', cc_data, 'MarkerSize', 15);
        h.sys(2,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit);
        
    ax(3,m) = subplot(3,n_plot,sI(3)) ; hold on ; box on
        %yline(1, '--k')
        h.data(3,m,:) = plot(IOFv, ALL.BodyFixed.FRF_data.C_adaptive(vI).error, ...
            '.-', 'Color', cc_data, 'MarkerSize', 15);
        h.sys(3,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit);
        ax(3).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(1:end,:), 'XLim', [0.1 20])
set(ax(1:end,:), 'XScale', 'log')
set(ax(1:end,:), 'XTick', [0.1 1 10])
set(ax(1,1), 'YLim', [0 3])
set(ax(2,:), 'YLim', [-100 50])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 6)
% set(h.sys_freq, 'MarkerSize', 6, 'LineWidth', 0.75)
% set(h.sys_freq, 'Marker', 'none')

set(ax(1:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Tracking error')

XLabelHC = get(ax(3,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

%% Closed-loop predicition
clc
MODEL.BodyFixed.H_predicition = MODEL.HeadFree.H.head_no_body;
for v = 1:3                                                
    ALL.BodyFixed.FRF_data.H_predicition(v).complex = ALL.HeadFree.FRF_data.err2head.grand_mean(v).complex ./ ...
                                                    (1 + ALL.HeadFree.FRF_data.err2head.grand_mean(v).complex); 
                                                
	ALL.BodyFixed.FRF_data.H_predicition(v).gain = abs(ALL.BodyFixed.FRF_data.H_predicition(v).complex);
    ALL.BodyFixed.FRF_data.H_predicition(v).phase = rad2deg( angle(ALL.BodyFixed.FRF_data.H_predicition(v).complex) );
    ALL.BodyFixed.FRF_data.H_predicition(v).error = abs( (1 + 0*1i) - ALL.BodyFixed.FRF_data.H_predicition(v).complex);
end

vI = 2;
D = {ALL.BodyFixed.FRF_data.H_predicition(vI).gain, ALL.BodyFixed.FRF_data.H_predicition(vI).phase};
M = plotFit(D, IOFv, MODEL.BodyFixed.H_predicition, 0:0.02:20, false);

clear ax h

cc_fit = [0 0.8 0.2];
cc_data = 'k';
n_plot = 1;

fig = figure (4); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [n_plot*(3/2) (5/4)*2.5];
movegui(fig, 'center')
ax = gobjects(3,n_plot);
n = 1;
for m = 1:n_plot
    sI = m + (0:2)*n_plot;        
    ax(1,m) = subplot(3,n_plot,sI(1)); cla ; hold on ; box on  
        h.sys(1,m) = plot(M(n).fv_bode, M(n).gain, 'Color', cc_fit);
        h.data(1,m,:) = plot(IOFv, D{1}, '.-', 'Color', cc_data, 'MarkerSize', 15);
        
    ax(2,m) = subplot(3,n_plot,sI(2)); cla ; hold on ; box on       
        h.sys(2,m) = plot(M(n).fv_bode, M(n).phase, 'Color', cc_fit);
        h.data(2,m,:) = plot(IOFv, D{2}, '.-', 'Color', cc_data, 'MarkerSize', 15);
        
    ax(3,m) = subplot(3,n_plot,sI(3)) ; hold on ; box on
        h.sys(3,m) = plot(M(n).fv_bode, M(n).error, 'Color', cc_fit);
        h.data(3,m,:) = plot(IOFv, ALL.BodyFixed.FRF_data.H_predicition(vI).error, ...
            '.-', 'Color', cc_data, 'MarkerSize', 15);
        ax(3).YLim(1) = 0;
end
set(ax, 'Color', 'none', 'LineWidth', 0.5)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax(1:end,:), 'XLim', [0.1 20])
set(ax(1:end,:), 'XScale', 'log')
set(ax(1:end,:), 'XTick', [0.1 1 10])
set(ax(1,1), 'YLim', [0 1])
set(ax(2,:), 'YLim', [-250 150])
linkaxes(ax(2:end,:), 'x')
set(h.data, 'MarkerSize', 6)

% set(h.sys_freq, 'MarkerSize', 6, 'LineWidth', 0.75)
% set(h.sys_freq, 'Marker', 'none')

set(ax(1:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

XLabelHC = get(ax(3,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

end