function [sys, data, h] = plant_power(P1, P2, showplot, fv)
%% plant_power: takes in transfer functions and compares there gains
%
%  INPUTS:
%       P1          : plant #1
%       P2          : plant #2
%
%  OUTPUTS:
%       sys         : structure containing model transfer functions
%      	data       	: bode plot and power data
%       h           : graphics handles
%

if nargin < 4
    fv = (0:0.05:30)';
    if nargin < 3
        showplot = true;
    end
end

% Plants
sys.P1 = P1; % plant #1
sys.P2 = P2; % plant #2

% Plant bode plot computations   
plot_tf = string(fieldnames(sys));    
n_tf = length(plot_tf);

data = [];
data.fv = fv;
data.win = data.fv * 2*pi;
for n = 1:n_tf
    [gain,phase,~] = bode(sys.(plot_tf(n)), data.win);
    if mean(phase > 200)
        phase = phase -360;
    end
    data.gain.(plot_tf(n)) = squeeze(gain);
    data.phase.(plot_tf(n)) = squeeze(phase);
end
data.power.ratio_1 = (data.gain.P2 - data.gain.P1) ./ data.gain.P1;

if showplot
    fig = figure (13);
    set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
    fig.Position(3:4) = [3 4];
    set(fig, 'Visible', 'on')
    clear ax
    
    ax(1) = subplot(2,1,1); hold on ; cla
        h.gain(1) = plot(data.fv, data.gain.P1, 'Color', 'r', 'LineWidth', 2);
        h.gain(2) = plot(data.fv, data.gain.P2, 'Color', 'b', 'LineWidth', 2);
        ylabel('Gain')

    ax(2) = subplot(2,1,2); hold on ; cla
        h.power(1) = plot(data.fv, 1 * data.power.ratio_1, 'Color', 'm', 'LineWidth', 2);
        xlabel('Frequency (hz')
        ylabel('Power factor')
        %ax(2).YLim(1) = -0.05;
    
    set(ax, 'Color', 'none', 'LineWidth', 1)
    set(ax, 'XScale', 'log', 'XLim', [0.1 data.fv(end)])
else
    h = [];    
end

end