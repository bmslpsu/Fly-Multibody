function [fcut] = compute_cutoff(fv, mag, phase, clss, showplot)
%% compute_cutoff: get approximate cutoff frequncies of higher order closed-loop responses
%
%  INPUTS:
%       fv       	: frequency vector
%       mag       	: magnitude data
%       phase     	: phase data
%       clss        : filter class "low" or "high"
%
%  OUTPUTS:
%       fcut      	: approximated frequency [hz]
%

if nargin < 5
    showplot = false;
end

switch clss
    case 'low'
        mag_cut = 0.7071*mag(1);
    case 'high'
        mag_cut = 0.7071*mag(end);
    otherwise
        error('must be "low" or "high"')
end

[~,fcutI] = min(abs(mag - mag_cut));
fcut = fv(fcutI);

if showplot
    fig = figure;
    set(fig, 'Color', 'w')
    ax(1) = subplot(2,1,1); cla ; hold on ; title(fcut)
        plot(fv, mag, 'b', 'LineWidth', 1.5)
        xline(fcut, '--r')
        yline(mag_cut, '--r')
    ax(2) = subplot(2,1,2); cla ; hold on
        plot(fv, phase, 'b', 'LineWidth', 1.5)

    set(ax, 'LineWidth', 1.5, 'Color', 'none')
    set(ax, 'XScale', 'log')
end

end

