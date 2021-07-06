function [sys, data, h] = plotFit(cmplx_rsp, freq, models, freq_range, showplot)
% plotFit: plots frequency response data and the response for the 
%          correspdonding fit tranfer function (s)
%
%  INPUTS:
%     cmplx_rsp : complex response, frequencies along rows & repetitions along columns
%          freq : frequency vector (hz)
%   	  model : transfer function or state-space models
%  	 freq_range : frequency range for simulation of models [freq_low freq high] (hz)
%      showplot : show figure boolean
%
%  OUTPUTS:
%      model : structure containing model parameters and fit TF model
%          h : structure containing figure grpahics objects
%

% Defaults
if nargin < 5
    showplot = true;
    if nargin < 4
        freq_range = [];
        if nargin < 3
            models = {};
        end
    end
end

if ~iscell(models)
    models = {models};
end

[n_freq, n_rep] = size(cmplx_rsp{1});

% System data
data.freq = freq;
data.freq_all = repmat(data.freq, [1 n_rep]);
data.input = cmplx_rsp;
data.gain = cmplx_rsp{1};
data.phase = cmplx_rsp{2};
% data.complex = data.gain .* exp(1i*data.phase);
data.complex = data.gain .* (cosd(data.phase) + 1i*sind(data.phase));
data.error = abs( (1 + 0*1i) - data.complex );
data.real = real(data.complex);
data.imag = imag(data.complex);

% Get frequency respnse for each model
n_model = length(models);
sys = [];
for n = 1:n_model
    M = models{n}; % current model
    sys(n).models = M;
    try
        sys(n).fitpercent.tf = M.Report.Fit.FitPercent;
    catch
        sys(n).fitpercent.tf = nan;
    end
    
    % Get the real & imaginary trajectory and gain & phase
    [R, I, wrad] = nyquist(M, 2*pi*freq_range);
    sys(n).real = squeeze(R);
    sys(n).imag = squeeze(I);
    sys(n).fv_nyquist = wrad ./ (2*pi);
    sys(n).error = abs( (1 + 0*1i) - complex(sys(n).real, sys(n).imag));
    
    [G, P, wrad] = bode(M, 2*pi*freq_range);
    sys(n).gain = squeeze(G);
    sys(n).phase = squeeze(P);
    sys(n).fv_bode = wrad / (2*pi);

    % Get the real & imaginary trajectory and gain & phase at input frequencies
    [R, I] = nyquist(M, 2*pi*data.freq); 
    sys(n).real_freq = squeeze(R);
    sys(n).imag_freq = squeeze(I);
    sys(n).error_freq = abs( (1 + 0*1i) - complex(sys(n).real_freq, sys(n).imag_freq));

  	[G, P] = bode(M, 2*pi*data.freq);
    sys(n).gain_freq = squeeze(G);
    sys(n).phase_freq = squeeze(P);
    
    phase_check = median(sys(n).phase_freq) - median(data.phase,'all');
    if abs(phase_check) > 180
        sys(n).phase_freq = sys(n).phase_freq - sign(phase_check)*360;
        sys(n).phase = sys(n).phase - sign(phase_check)*360;
    end

    % Compute goodness of fit metrics
    gain_freq_all = repmat(sys(n).gain_freq, [n_rep 1]);
    phase_freq_all = repmat(sys(n).phase_freq, [n_rep 1]);
    gain_nmse = goodnessOfFit(data.gain(:), gain_freq_all, 'NMSE');
    phase_nmse = goodnessOfFit(data.phase(:), phase_freq_all, 'NMSE');
    
    %sys(n).fitpercent.gain = 1 - gain_nmse;
    %sys(n).fitpercent.phase = 1 - phase_nmse;
    sys(n).fitpercent.gain = gain_nmse;
    sys(n).fitpercent.phase = phase_nmse;
    sys(n).fitpercent.combined = mean([sys(n).fitpercent.gain sys(n).fitpercent.phase]);
end

if showplot
    if isa(showplot,'double')
        cc = repmat([0.5 0.5 0.5], [n_freq, 1]);
        cc_fit = lines(n_model);
    else
        cc = hsv(n_freq);
        cc_fit = gray(n_model + 1);
        cc_fit = cc_fit(1:n_model,:);
    end
    fig = figure; clf
    set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
    fig.Position(3:4) = [4 2.5*4];
    movegui(fig, 'center')
    ax(1) = subplot(4,1,1); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k');
        yline(0, '--k');
        %plot(1, 0, 'g.', 'MarkerSize', 25)
        for n = 1:n_model
            h.data(1,n,:) = gscatter(data.real(:), data.imag(:), data.freq_all(:), ...
                cc, '.', 10, false);
            h.sys(1,n) = plot(sys(n).real, sys(n).imag, 'Color', cc_fit(n,:));
            h.sys_freq(1,n) = plot(sys(n).real_freq, sys(n).imag_freq, ...
                '.', 'Color', cc_fit(n,:), 'MarkerSize', 15);
            plot(sys(n).real(1), sys(n).imag(1), '.', 'Color', 'm', 'MarkerSize', 15)
        end
        xlabel('Real')
        ylabel('Imaginary')
    ax(2) = subplot(4,1,2); cla ; hold on ; box on
        for n = 1:n_model
            h.data(2,n,:) = gscatter(data.freq_all(:), data.gain(:), data.freq_all(:), ...
                cc, '.', 10, false);
            h.sys(2,n) = plot(sys(n).fv_bode, sys(n).gain, 'Color', cc_fit(n,:));
            h.sys_freq(2,n) = plot(data.freq, sys(n).gain_freq, ...
                '.', 'Color', cc_fit(n,:), 'MarkerSize', 15);
        end
        ylabel('Gain')
    ax(3) = subplot(4,1,3); cla ; hold on ; box on
        for n = 1:n_model
            h.data(3,n,:) = gscatter(data.freq_all(:), data.phase(:), data.freq_all(:), ...
                cc, '.', 10, false);
            h.sys(3,n) = plot(sys(n).fv_bode, sys(n).phase, 'Color', cc_fit(n,:));
            h.sys_freq(3,n) = plot(data.freq, sys(n).phase_freq, ...
                '.', 'Color', cc_fit(n,:), 'MarkerSize', 15);
        end
        ylabel('Phase')
        
    ax(4) = subplot(4,1,4); cla ; hold on ; box on
        yline(1, '--k');
        for n = 1:n_model
            h.data(4,n,:) = gscatter(data.freq_all(:), data.error(:), data.freq_all(:), ...
                cc, '.', 10, false);
            h.sys(4,n) = plot(sys(n).fv_bode, sys(n).error, 'Color', cc_fit(n,:));
            h.sys_freq(4,n) = plot(data.freq, sys(n).error_freq, ...
                '.', 'Color', cc_fit(n,:), 'MarkerSize', 15);
        end
        ax(4).YLim(1) = 0;
        xlabel('Frequency (Hz)')
        ylabel('Tracking error')

    set(ax, 'Color', 'none', 'LineWidth', 1)
    set(ax, 'XGrid', 'on', 'YGrid', 'on')
    set(ax(2:end), 'XLim', [0.1 20])
    set(ax(2:end), 'XScale', 'log')
    linkaxes(ax(2:end), 'x')

    h.fig = fig;
    h.ax = ax;
else
    h = [];
end

end