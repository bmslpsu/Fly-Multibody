function [] = SOS_power()
%% SOS_power:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

load('E:\DATA\Magno_Data\Multibody\sym\Symbolic_expressions.mat')

%% Set initialization parameters
clearvars -except PATH FILE ALL MODEL equations
clc

P_body = MODEL.HeadFree.P.body;
C_body = MODEL.HeadFree.C.body;

P_head = MODEL.HeadFree.P.head;
C_head = MODEL.HeadFree.C.head;

G_body = P_body * C_body;
G_head = P_head * C_head;

H_body = minreal( G_body / (1 + G_body + G_head) );
H_head = minreal( G_head / (1 + G_body + G_head) );
H_gaze = minreal(H_body + H_head);

W_body = minreal( C_body / (1 + G_body + G_head) );
W_head = minreal( C_head / (1 + G_body + G_head) );
W_gaze = minreal(W_body + W_head);

% h = bodeplot(H_body, H_head, H_gaze, MODEL.HeadFree.H.body, MODEL.HeadFree.H.head, MODEL.HeadFree.H.gaze);
% setoptions(h,'FreqUnits','Hz', 'MagUnits', 'abs')
% xlim([0.1 20])

syms s H1 P1 C1
C_body_all = subs(equations.single.controllers.C1, [H1 P1], [tf2sym(H_gaze) tf2sym(P_body)]);
C_body_all = minreal_sym(collect(simplify(expand(C_body_all)),s),s);

H_body_all = subs(equations.single.closedloop.H1, [C1 P1], [C_body_all tf2sym(P_body)]);
H_body_all = minreal_sym(collect(simplify(expand(H_body_all)),s),s);
H_body_all = minreal(sym2tf(H_body_all));

C_body_all = minreal(sym2tf(C_body_all));
G_body_all = P_body * C_body_all;

W_body_all = minreal( C_body_all / (1 + G_body_all) );

%%
plot_tf = {P_body, P_head, C_body, C_head, C_body_all, ...
            H_body, H_head, H_gaze, H_body_all, W_body, W_head, W_body_all};
        
n_plot = length(plot_tf);
cc = [1 0 0 ; 0 0 1 ; 0.3 0.2 0.6 ; 0 1 0.7];
fv = (0:0.001:20)';
win = fv * 2*pi;
n_freq = length(fv);
Mag = nan(n_freq, n_plot);
Phs = nan(n_freq, n_plot);
for n = 1:n_plot
    [mag,phase,~] = bode(plot_tf{n}, win);
    Mag(:,n) = squeeze(mag);
    Phs(:,n) = squeeze(phase);
end

fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'simulation')
fig.Position(3:4) = [9 4.5];
set(fig, 'Visible', 'on')
clear ax h
ax = gobjects(2,4);
ax(1,1) = subplot(2,4,1); cla ; hold on ; title('Plant')
    h.plant(:,1) = plot(fv, Mag(:,1:2));
    ylabel('gain (°/°)')
    
ax(2,1) = subplot(2,4,5); cla ; hold on
    yline(0, '--k');
    h.plant(:,2) = plot(fv, Phs(:,1:2));
    ylabel('phase (°)')
    xlabel('frequency (hz)')
    
ax(1,2) = subplot(2,4,2); cla ; hold on ; title('Controller')
    h.cntrl(:,1) = plot(fv, Mag(:,3:5));
    
ax(2,2) = subplot(2,4,6); cla ; hold on
    yline(0, '--k');
    h.cntrl(:,2) = plot(fv, Phs(:,3:5));
    xlabel('frequency (hz)')
    
ax(1,3) = subplot(2,4,3); cla ; hold on ; title('Closed-loop')
    h.cl(:,1) = plot(fv, Mag(:,6:9));
    
ax(2,3) = subplot(2,4,7); cla ; hold on
    yline(0, '--k');
    h.cl(:,2) = plot(fv, Phs(:,6:9));
    xlabel('frequency (hz)')
    
ax(1,4) = subplot(2,4,4); cla ; hold on ; title('Control input')
    h.u(:,1) = plot(fv, Mag(:,10:12));
	h.u_sum(:,1) = plot(fv, Mag(:,10) + Mag(:,11));
    
ax(2,4) = subplot(2,4,8); cla ; hold on
    yline(0, '--k');
    h.u(:,2) = plot(fv, Phs(:,10:12));
    xlabel('frequency (hz)')

linkaxes(ax, 'x')
linkaxes(ax(2,:), 'y')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XScale', 'log')
set(ax, 'XLim', [0.1 20], 'XTick', [0.1 1 10])
% set(ax(1,:), 'XColor', 'none')
set(h.cl(end,:), 'LineStyle', '--')

structfun(@(x) set(x, 'LineWidth', 1.5), h)
set(h.plant(:,1),{'Color'}, num2cell(cc(1:2,:),2))
set(h.plant(:,2),{'Color'}, num2cell(cc(1:2,:),2))
set(h.cntrl(:,1),{'Color'}, num2cell(cc([1:2,4],:),2))
set(h.cntrl(:,2),{'Color'}, num2cell(cc([1:2,4],:),2))
set(h.cl(:,1),{'Color'}, num2cell(cc(1:4,:),2))
set(h.cl(:,2),{'Color'}, num2cell(cc(1:4,:),2))
set(h.u(:,1),{'Color'}, num2cell(cc([1:2,4],:),2))
set(h.u(:,2),{'Color'}, num2cell(cc([1:2,4],:),2))
set(h.u_sum(:,1),{'Color'}, num2cell(cc(3,:),2))

end