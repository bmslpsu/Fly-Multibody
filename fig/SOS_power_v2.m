function [] = SOS_power_v2()
%% SOS_power:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

load('E:\DATA\Magno_Data\Multibody\sym\Symbolic_expressions.mat')

%% 
clearvars -except PATH FILE ALL MODEL equations
clc

% sys.P_body = MODEL.HeadFree.P.body;
% sys.C_body = MODEL.HeadFree.C.body;
% 
% sys.P_head = MODEL.HeadFree.P.head;
% sys.C_head = MODEL.HeadFree.C.head;

sys.P_body = MODEL.HeadFree.P_norm.body;
sys.C_body = MODEL.HeadFree.C_norm.body;

sys.P_head = MODEL.HeadFree.P_norm.head;
sys.C_head = MODEL.HeadFree.C_norm.head;

sys.G_body = sys.P_body * sys.C_body;
sys.G_head = sys.P_head * sys.C_head;
sys.G_body.IODelay = MODEL.HeadFree.G.body.IODelay;
sys.G_head.IODelay = MODEL.HeadFree.G.head.IODelay;
sys.G_body = pade(sys.G_body);
sys.G_head = pade(sys.G_head);

sys.H_body = minreal( sys.G_body / (1 + sys.G_body + sys.G_head) );
sys.H_head = minreal( sys.G_head / (1 + sys.G_body + sys.G_head) );
sys.H_gaze = minreal( sys.H_body + sys.H_head );

sys.W_body = minreal( sys.C_body / (1 + sys.G_body + sys.G_head) );
sys.W_head = minreal( sys.C_head / (1 + sys.G_body + sys.G_head) );
sys.W_gaze = minreal( sys.W_body + sys.W_head );

% Simulated body response
syms s H1 P1 C1
sys.C_body_all = subs(equations.single.controllers.C1, [H1 P1], [tf2sym(sys.H_head) tf2sym(sys.P_body)]);
sys.C_body_all = minreal_sym(collect(simplify(expand(sys.C_body_all)),s),s);

sys.H_body_all = subs(equations.single.closedloop.H1, [C1 P1], [sys.C_body_all tf2sym(sys.P_body)]);
sys.H_body_all = minreal_sym(collect(simplify(expand(sys.H_body_all)),s),s);
sys.H_body_all = minreal(sym2tf(sys.H_body_all));

sys.C_body_all = minreal( sym2tf(sys.C_body_all) );
sys.G_body_all = sys.P_body * sys.C_body_all;

sys.W_body_all = minreal( sys.C_body_all / (1 + sys.G_body_all) );

% Simulated head response
sys.C_head_all = subs(equations.single.controllers.C1, [H1 P1], [tf2sym(sys.H_gaze) tf2sym(sys.P_head)]);
sys.C_head_all = minreal_sym(collect(simplify(expand(sys.C_head_all)),s),s);

sys.H_head_all = subs(equations.single.closedloop.H1, [C1 P1], [sys.C_head_all tf2sym(sys.P_head)]);
sys.H_head_all = minreal_sym(collect(simplify(expand(sys.H_head_all)),s),s);
sys.H_head_all = minreal(sym2tf(sys.H_head_all));

sys.C_head_all = minreal( sym2tf(sys.C_head_all) );
sys.G_head_all = sys.P_head * sys.C_head_all;

sys.W_head_all = minreal( sys.C_body_all / (1 + sys.G_body_all) );

%% Plant, Conctroller, CLosed-loop, Control-input
clc
plot_tf = ["P_body", "P_head", "C_body", "C_head", "C_body_all", "C_head_all"...
           "H_body", "H_head", "H_gaze", "H_body_all", "H_head_all", ...
           "W_body", "W_head", "W_gaze", "W_body_all", "W_head_all"];
        
n_plot = length(plot_tf);

data = [];
data.fv = (0:0.001:20)';
data.iofv = ALL.HeadFree.FUNC{2, 1}.All.Freq;
data.win = data.fv * 2*pi;
data.iowin = data.iofv * 2*pi;
for n = 1:n_plot
    [gain,phase,~] = bode(sys.(plot_tf(n)), data.win);
    if mean(phase > 200)
        phase = phase -360;
    end
    data.gain.(plot_tf(n)) = squeeze(gain);
    data.phase.(plot_tf(n)) = squeeze(phase);
    
    [gain,phase,~] = bode(sys.(plot_tf(n)), data.iowin);
    if mean(phase > 200)
        phase = phase -360;
    end
    data.gain_io.(plot_tf(n)) = squeeze(gain);
    data.phase_io.(plot_tf(n)) = squeeze(phase);
    
    name = char(plot_tf(n));
    if strcmp(name(1),'W')
        data.power.(plot_tf(n)) = data.fv .* data.gain.(plot_tf(n));
        data.power_io.(plot_tf(n)) = data.iofv .* data.gain_io.(plot_tf(n));
    end
end
data.gain.Csum = data.gain.C_body + data.gain.C_head;
data.gain.Wsum = data.gain.W_body + data.gain.W_head;
data.power.Wsum = data.power.W_body + data.power.W_head;
data.power.ratio = data.power.W_head ./ data.power.W_body_all;

data.gain_io.Csum = data.gain_io.C_body + data.gain_io.C_head;
data.gain_io.Wsum = data.gain_io.W_body + data.gain_io.W_head;
data.power_io.Wsum = data.power_io.W_body + data.power_io.W_head;
data.power_io.ratio = data.power_io.W_head ./ data.power_io.W_body_all;

p = data.power.ratio;
mp = max(p);
p(data.fv < mp) = nan;
[~,fI] = min(abs(p - 1));
data.power.freq_improve = data.fv(fI);

%%
cc = [];
cc.body = [0.9 0 0];
cc.head = [0 0.4 1];
cc.gaze = [0.5 0.3 1];
cc.gaze_comb = [0 1 0.7];

fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [9 4.5];
set(fig, 'Visible', 'on')
clear ax h
ax = gobjects(2,4);
ax(1,1) = subplot(2,4,1); cla ; hold on ; title('Plant') ; ylabel('gain')
    h.P(1,1) = plot(data.fv, data.gain.P_body, 'Color', cc.body);
    h.P(1,2) = plot(data.fv, data.gain.P_head, 'Color', cc.head);
    
   	h.P(2,1) = plot(data.iofv, data.gain_io.P_body, 'Color', cc.body);
    h.P(2,2) = plot(data.iofv, data.gain_io.P_head, 'Color', cc.head);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
        
ax(2,1) = subplot(2,4,5); cla ; hold on ; ylabel('phase (°)')
    yline(0, '--k')
    h.P(1,3) = plot(data.fv, data.phase.P_body, 'Color', cc.body);
    h.P(1,4) = plot(data.fv, data.phase.P_head, 'Color', cc.head);
    
  	h.P(2,3) = plot(data.iofv, data.phase_io.P_body, 'Color', cc.body);
    h.P(2,4) = plot(data.iofv, data.phase_io.P_head, 'Color', cc.head);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
        
ax(1,2) = subplot(2,4,2); cla ; hold on ; title('Controller')
    h.C(1,1) = plot(data.fv, data.gain.C_body, 'Color', cc.body);
    h.C(1,2) = plot(data.fv, data.gain.C_head, 'Color', cc.head);
    h.C(1,3) = plot(data.fv, data.gain.Csum, 'Color', cc.gaze);
    h.C(1,4) = plot(data.fv, data.gain.C_body_all, 'Color', cc.gaze_comb);
    
    h.C(2,1) = plot(data.iofv, data.gain_io.C_body, 'Color', cc.body);
    h.C(2,2) = plot(data.iofv, data.gain_io.C_head, 'Color', cc.head);
    h.C(2,3) = plot(data.iofv, data.gain_io.Csum, 'Color', cc.gaze);
    h.C(2,4) = plot(data.iofv, data.gain_io.C_body_all, 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
    
ax(2,2) = subplot(2,4,6); cla ; hold on
    yline(0, '--k')
    h.C(1,5) = plot(data.fv, data.phase.C_body, 'Color', cc.body);
    h.C(1,6) = plot(data.fv, data.phase.C_head, 'Color', cc.head);
    h.C(1,7) = plot(data.fv, data.phase.C_body_all, 'Color', cc.gaze_comb);
    
    h.C(2,5) = plot(data.iofv, data.phase_io.C_body, 'Color', cc.body);
    h.C(2,6) = plot(data.iofv, data.phase_io.C_head, 'Color', cc.head);
    h.C(2,7) = plot(data.iofv, data.phase_io.C_body_all, 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
    
ax(1,3) = subplot(2,4,3); cla ; hold on ; title('Closed-loop')
    h.H(1,1) = plot(data.fv, data.gain.H_body, 'Color', cc.body);
    h.H(1,2) = plot(data.fv, data.gain.H_head, 'Color', cc.head);
    h.H(1,3) = plot(data.fv, data.gain.H_gaze, 'Color', cc.gaze);
    h.H(1,4) = plot(data.fv, data.gain.H_body_all, '--', 'Color', cc.gaze_comb);
    
    h.H(2,1) = plot(data.iofv, data.gain_io.H_body, 'Color', cc.body);
    h.H(2,2) = plot(data.iofv, data.gain_io.H_head, 'Color', cc.head);
    h.H(2,3) = plot(data.iofv, data.gain_io.H_gaze, 'Color', cc.gaze);
    h.H(2,4) = plot(data.iofv, data.gain_io.H_body_all, '--', 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
    
ax(2,3) = subplot(2,4,7); cla ; hold on
    yline(0, '--k')
    h.H(1,5) = plot(data.fv, data.phase.H_body, 'Color', cc.body);
    h.H(1,6) = plot(data.fv, data.phase.H_head, 'Color', cc.head);
    h.H(1,7) = plot(data.fv, data.phase.H_gaze, 'Color', cc.gaze);
    h.H(1,8) = plot(data.fv, data.phase.H_body_all, '--', 'Color', cc.gaze_comb);
    
    h.H(2,5) = plot(data.iofv, data.phase_io.H_body, 'Color', cc.body);
    h.H(2,6) = plot(data.iofv, data.phase_io.H_head, 'Color', cc.head);
    h.H(2,7) = plot(data.iofv, data.phase_io.H_gaze, 'Color', cc.gaze);
    h.H(2,8) = plot(data.iofv, data.phase_io.H_body_all, '--', 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
    
ax(1,4) = subplot(2,4,4); cla ; hold on ; title('Control input')
    h.W(1,1) = plot(data.fv, data.gain.W_body, 'Color', cc.body);
    h.W(1,2) = plot(data.fv, data.gain.W_head, 'Color', cc.head);
    h.W(1,3) = plot(data.fv, data.gain.Wsum, 'Color', cc.gaze);
    h.W(1,4) = plot(data.fv, data.gain.W_body_all, 'Color', cc.gaze_comb);
    
  	h.W(2,1) = plot(data.iofv, data.gain_io.W_body, 'Color', cc.body);
    h.W(2,2) = plot(data.iofv, data.gain_io.W_head, 'Color', cc.head);
    h.W(2,3) = plot(data.iofv, data.gain_io.Wsum, 'Color', cc.gaze);
    h.W(2,4) = plot(data.iofv, data.gain_io.W_body_all, 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
    
ax(2,4) = subplot(2,4,8); cla ; hold on
    yline(0, '--k')
    h.W(1,5) = plot(data.fv, data.phase.W_body, 'Color', cc.body);
    h.W(1,6) = plot(data.fv, data.phase.W_head, 'Color', cc.head);
    h.W(1,7) = plot(data.fv, data.phase.W_body_all, 'Color', cc.gaze_comb);
    
    h.W(2,5) = plot(data.iofv, data.phase_io.W_body, 'Color', cc.body);
    h.W(2,6) = plot(data.iofv, data.phase_io.W_head, 'Color', cc.head);
    h.W(2,7) = plot(data.iofv, data.phase_io.W_body_all, 'Color', cc.gaze_comb);
    xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])

linkaxes(ax, 'x')
linkaxes(ax(2,:), 'y')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XScale', 'log')
set(ax, 'XLim', [0.1 20], 'XTick', [0.1 1 10])
set(ax(1,:), 'XtickLabels', [])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', {'frequency (hz)'})

structfun(@(x) set(x, 'LineWidth', 1.5), h)
structfun(@(x) set(x(2,:), 'Marker', 'none', 'MarkerSize', 12, 'LineStyle', 'none'), h)

%% 
fig = figure (2);
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [3 2];
set(fig, 'Visible', 'on')
clear ax h

% cc = [0.8 0.3 0.5];
cc = [0 0.8 0.3];
ax = subplot(1,1,1); hold on ; cla
yline(1, 'k--')
plot(data.fv, data.power.ratio, 'Color', cc, 'LineWidth', 2)
plot(data.iofv, data.power_io.ratio, '.', 'Color', cc, 'MarkerSize', 12)
xline(data.power.freq_improve, '--', 'Color', [0.5 0.5 0.5])
ylim([0 1.5])
xlabel('Frequency (hz')
ylabel('Power factor')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax, 'XScale', 'log')
xlim([0.1 20])

end