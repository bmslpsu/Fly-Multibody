function [sys, data, h] = multi_body_system(P1, P2, C1, C2, delay1, delay2, colors)
% multi_body_system: takes in parallel plant and controller transfer functions
% (& optional correspondiong delays) and computes the closed-loop (H) & control input (W)
% transfer functions.
%
%  INPUTS:
%       P1          : plant #1
%       P2          : plant #2
%       C1          : controller #1
%       C2          : controller #2
%       delay1      : optional iodelay for P1*C1
%       delay2      : optional iodelay for P2*C2
%       showplot    : show figure boolean
%
%  OUTPUTS:
%       sys         : structure containing model transfer functions
%      	data       	: bode plot and power data
%       h           : graphics handles
%

if nargin < 7
    colors = [];
    if nargin < 6 
        delay2 = 0;
        if nargin < 5
            delay1 = 0;
        end
    end
end

% 1st open-loop system
sys.P1 = P1; % plant #1
sys.C1 = C1; % controller #1
sys.G1 = sys.P1 * sys.C1; % open-loop system #1
sys.G1.IODelay = delay1; % set delay
sys.G1 = pade(sys.G1); % approximate as linear system if delay is nonzero

% 2nd open-loop system
sys.P2 = P2; % plant #2
sys.C2 = C2; % controller #2
sys.G2 = sys.P2 * sys.C2; % open-loop system #2
sys.G2.IODelay = delay2; % set delay
sys.G2 = pade(sys.G2); % approximate as linear system if delay is nonzero

% Closed-loop systems and combined system
sys.H1 = minreal( sys.G1 / (1 + sys.G1 + sys.G2) ); % closed-loop system #1
sys.H2 = minreal( sys.G2 / (1 + sys.G1 + sys.G2) ); % closed-loop system #2
sys.H12 = minreal( sys.H1 + sys.H2 ); % full closed-loop system

% Control input transfer functions
sys.W1 = minreal( sys.C1 / (1 + sys.G1 + sys.G2) ); % control-input system #1
sys.W2 = minreal( sys.C2 / (1 + sys.G1 + sys.G2) ); % control-input system #2

%% Plant, Controller, Closed-loop, Control-input bode plot computations   
plot_tf = string(fieldnames(sys));    
n_tf = length(plot_tf);

data = [];
data.fv = (0:0.001:100)';
data.win = data.fv * 2*pi;
for n = 1:n_tf
    [gain,phase,~] = bode(sys.(plot_tf(n)), data.win);
    if mean(phase > 200)
        phase = phase -360;
    end
    data.gain.(plot_tf(n)) = squeeze(gain);
    data.phase.(plot_tf(n)) = squeeze(phase);
    
    name = char(plot_tf(n));
    if strcmp(name(1),'W') || strcmp(name(1),'C')
        data.power.(plot_tf(n)) = data.fv .* data.gain.(plot_tf(n));
    end
end
data.gain.Csum = data.gain.C1 + data.gain.C2;
data.gain.Wsum = data.gain.W1 + data.gain.W2;
data.power.Wsum = data.power.W1 + data.power.W2;
% data.power.ratio_1 = data.power.W2 ./ data.power.W1_switch;
% data.power.ratio_2 = data.power.W1 ./ data.power.W2_switch;
data.power.ratio_1 = (data.power.W1_switch - data.power.W2) ./ data.power.W2;
data.power.ratio_2 = (data.power.W2_switch - data.power.W1) ./ data.power.W1;

% data.power.ratio_1 = (data.power.C1_switch - data.power.C2) ./ data.power.C2;
% data.power.ratio_2 = (data.power.C2_switch - data.power.C1) ./ data.power.C1;

%% Plots
if ~isempty(colors)
    if islogical(colors)
        if colors 
            showplot = colors;
            
            disp('Using default colors')
            colors = {};
          	colors{1} = [0.9 0 0];
            colors{2} = [0 0.4 1];
            colors{3} = [0.5 0.3 1];
            colors{4} = [0.9 0.4 0.17];
            colors{5} = [0 1 0.7];
        else
            showplot = false;
            h = [];
        end
    else
        showplot = true;
        disp('Setting colors manually')
    end
end

if showplot
cc.one = colors{1};
cc.two = colors{2};
cc.comb = colors{3};
cc.switch1 = colors{4};
cc.switch2 = colors{5};

fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [9 4.5];
set(fig, 'Visible', 'on')
clear ax h
ax = gobjects(2,4);
ax(1,1) = subplot(2,4,1); cla ; hold on ; title('Plant') ; ylabel('gain')
    h.P(1,1) = plot(data.fv, data.gain.P1, 'Color', cc.one);
    h.P(1,2) = plot(data.fv, data.gain.P2, 'Color', cc.two);
            
ax(2,1) = subplot(2,4,5); cla ; hold on ; ylabel('phase (°)')
    yline(0, '--k');
    h.P(1,3) = plot(data.fv, data.phase.P1, 'Color', cc.one);
    h.P(1,4) = plot(data.fv, data.phase.P2, 'Color', cc.two);
            
ax(1,2) = subplot(2,4,2); cla ; hold on ; title('Controller')
    h.C(1,1) = plot(data.fv, data.gain.C1, 'Color', cc.one);
    h.C(1,2) = plot(data.fv, data.gain.C2, 'Color', cc.two);
    h.C(1,3) = plot(data.fv, data.gain.Csum, 'Color', cc.comb);
    h.C(1,4) = plot(data.fv, data.gain.C1_switch, '--', 'Color', cc.switch1);
    h.C(1,5) = plot(data.fv, data.gain.C2_switch, '--', 'Color', cc.switch2);
        
ax(2,2) = subplot(2,4,6); cla ; hold on
    yline(0, '--k');
    h.C(1,6) = plot(data.fv, data.phase.C1, 'Color', cc.one);
    h.C(1,7) = plot(data.fv, data.phase.C2, 'Color', cc.two);
    h.C(1,8) = plot(data.fv, data.phase.C1_switch, '--', 'Color', cc.switch1);
    h.C(1,9) = plot(data.fv, data.phase.C2_switch, '--', 'Color', cc.switch2);
        
ax(1,3) = subplot(2,4,3); cla ; hold on ; title('Closed-loop')
    h.H(1,1) = plot(data.fv, data.gain.H1, 'Color', cc.one);
    h.H(1,2) = plot(data.fv, data.gain.H2, 'Color', cc.two);
    h.H(1,3) = plot(data.fv, data.gain.H12, 'Color', cc.comb);
    h.H(1,4) = plot(data.fv, data.gain.H1_switch, '--', 'Color', cc.switch1);
  	h.H(1,5) = plot(data.fv, data.gain.H2_switch, '--', 'Color', cc.switch2);
        
ax(2,3) = subplot(2,4,7); cla ; hold on
    yline(0, '--k');
    h.H(1,6) = plot(data.fv, data.phase.H1, 'Color', cc.one);
    h.H(1,7) = plot(data.fv, data.phase.H2, 'Color', cc.two);
    h.H(1,8) = plot(data.fv, data.phase.H12, 'Color', cc.comb);
    h.H(1,9) = plot(data.fv, data.phase.H1_switch, '--', 'Color', cc.switch1);
    h.H(1,10) = plot(data.fv, data.phase.H2_switch, '--', 'Color', cc.switch2);
        
ax(1,4) = subplot(2,4,4); cla ; hold on ; title('Control input')
    h.W(1,1) = plot(data.fv, data.gain.W1, 'Color', cc.one);
    h.W(1,2) = plot(data.fv, data.gain.W2, 'Color', cc.two);
    h.W(1,3) = plot(data.fv, data.gain.Wsum, 'Color', cc.comb);
    h.W(1,4) = plot(data.fv, data.gain.W1_switch, '--', 'Color', cc.switch1);
    h.W(1,5) = plot(data.fv, data.gain.W2_switch, '--', 'Color', cc.switch2);
        
ax(2,4) = subplot(2,4,8); cla ; hold on
    yline(0, '--k');
    h.W(1,6) = plot(data.fv, data.phase.W1, 'Color', cc.one);
    h.W(1,7) = plot(data.fv, data.phase.W2, 'Color', cc.two);
    h.W(1,8) = plot(data.fv, data.phase.W1_switch, '--', 'Color', cc.switch1);
    h.W(1,9) = plot(data.fv, data.phase.W2_switch, '--', 'Color', cc.switch2);

linkaxes(ax, 'x')
linkaxes(ax(2,:), 'y')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XScale', 'log')
set(ax, 'XLim', [0.1 20], 'XTick', [0.1 1 10])
set(ax(1,:), 'XtickLabels', [])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', {'frequency (hz)'})

structfun(@(x) set(x, 'LineWidth', 1.5), h)

%% Power
fig = figure (2);
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [3 2];
set(fig, 'Visible', 'on')
clear ax

% cc = [0.6 0.6 0.6];
ax = subplot(1,1,1); hold on ; cla
yline(1, 'k--');
h.power(1) = plot(data.fv, data.power.ratio_1, 'Color', cc.switch1, 'LineWidth', 2);
% h.power(2) = plot(data.fv, data.power.ratio_2, 'Color', cc.switch2, 'LineWidth', 2);
xlabel('Frequency (hz')
ylabel('Power factor')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax, 'XScale', 'log')
xlim([0.1 20])
ax.YLim(1) = -0.05;

end

end