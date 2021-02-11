function [] = SS_frf()
%% SS_frf:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','FUNC','GRAND','FLY','D','I','U','N')

%%
clc
clearvars -except FILE DATA ALL GRAND FLY FUNC D I U N root

pI = [1 2 3 4 5 6 7 8];
T = ["ref2body", "ref2head", "ref2gaze", "head2body", "ref2wing", "wing2body", "error2body", "error2head"];
shift_I = {5:8, 7:8, 5:8, 5:8, 5:8, 5:8, 5:8, 7:8};
phase_lim = {[0 0 -50 -220], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [nan nan]};

% pI = [1 2 3];
% T = ["ref2head", "ref2wing", "head2wing"];
% shift_I = {5:8, 5:8, 5:8, 5:8, 5:8, 5:8};
% phase_lim = {[0 0 0 -120], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0]};

% pI = [1 2 3 4];
% T = ["ref2body", "ref2head", "ref2gaze", "head2body"];
% shift_I = {5:8, 7:8, 5:8, 5:8};
% phase_lim = {[0 0 0 0], nan*[0 0 0 0], nan*[0 0 0 0], [0 0 0 0]};

n_plot = length(pI);
for n = 1:n_plot
    for v = 1:N.freq
        % Frequency components
        FRF_data.IOFv(v,1) = GRAND.fly_stats(v).mean.IOFv.mean;
        FRF_data.Fv = GRAND.fly_stats(v).mean.Fv.mean;
        
        % Input
        FRF_data.ref_mag(v) = GRAND.fly_stats(v).mean.refIOMag.mean(1);
        
        % FRF properties
        mag_all = squeeze(GRAND.fly_all(v).mean.IOMag(:,pI(n),:));
        mag_med = GRAND.fly_stats(v).mean.IOMag.mean(1,pI(n));
        mag_std = GRAND.fly_stats(v).mean.IOMag.std(1,pI(n));

        gain_all = squeeze(GRAND.fly_all(v).mean.IOGain(:,pI(n),:));
        gain_med = GRAND.fly_stats(v).mean.IOGain.mean(1,pI(n));
        gain_std = GRAND.fly_stats(v).mean.IOGain.std(1,pI(n));

        phase_all = rad2deg(squeeze(GRAND.fly_all(v).circ_mean.IOPhaseDiff(:,pI(n),:)));
        for r = 1:3
            shift_all = any(v == shift_I{n}, 2) & (phase_all > phase_lim{n}(v == shift_I{n}));
            phase_all(shift_all) = phase_all(shift_all) - 360;
        end
        
        phase_med = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_mean(:,pI(n)));
        phase_std = rad2deg(GRAND.fly_stats(v).circ_mean.IOPhaseDiff.circ_std(:,pI(n)));
        
        phase_med = mean(phase_all);
        phase_std = std(phase_all);
        for r = 1:2
            shift_all = any(v == shift_I{n}, 2) & (phase_med > phase_lim{n}(v == shift_I{n}));
            phase_med(shift_all) = phase_med(shift_all) - 360;
        end

        time_diff_all = 1000 * (phase_all ./360) .* (1 ./FRF_data.IOFv(v));
        time_diff_med = 1000 * (phase_med ./360) .* (1 ./FRF_data.IOFv(v));
        time_diff_std = 1000 * (phase_std ./360) .* (1 ./FRF_data.IOFv(v));

        cohr_all = squeeze(GRAND.fly_all(v).mean.Cohr(:,pI(n),:));
        cohr_med = GRAND.fly_stats(v).mean.Cohr.mean(:,pI(n));
        cohr_std = GRAND.fly_stats(v).mean.Cohr.std(:,pI(n));
        
      	IOcohr_all = squeeze(GRAND.fly_all(v).mean.IOCohr(:,pI(n),:));
        IOcohr_med = GRAND.fly_stats(v).mean.IOCohr.mean(:,pI(n));
        IOcohr_std = GRAND.fly_stats(v).mean.IOCohr.std(:,pI(n));
        
        % Fly means
        FRF_data.(T(n)).fly.mag(v,:) = mag_all;
        FRF_data.(T(n)).fly.gain(v,:) = gain_all;
        FRF_data.(T(n)).fly.phase(v,:) = phase_all;
        FRF_data.(T(n)).fly.time_diff(v,:) = time_diff_all;
        FRF_data.(T(n)).fly.coherence{v,1} = cohr_all;
        FRF_data.(T(n)).fly.IO_coherence(v,:) = IOcohr_all;
     	FRF_data.(T(n)).fly.complex(v,:) = gain_all .* (cosd(phase_all) + 1i*sind(phase_all));
        FRF_data.(T(n)).fly.error(v,:) = abs((1 + 0*1i) - FRF_data.(T(n)).fly.complex(v,:));
        
        % Grand means
        FRF_data.(T(n)).grand_mean.mag(v,1) = mag_med;
        FRF_data.(T(n)).grand_mean.gain(v,1) = gain_med;
        FRF_data.(T(n)).grand_mean.phase(v,1) = phase_med;
        FRF_data.(T(n)).grand_mean.time_diff(v,1) = time_diff_med;
        FRF_data.(T(n)).grand_mean.coherence(:,v) = cohr_med;
        FRF_data.(T(n)).grand_mean.IO_coherence(v,1) = IOcohr_med;
        FRF_data.(T(n)).grand_mean.complex(v,1) = mean(FRF_data.(T(n)).fly.complex(v,:));
        FRF_data.(T(n)).grand_mean.error(v,1) = mean(FRF_data.(T(n)).fly.error(v,:));
        
        % Grand STD's
        FRF_data.(T(n)).grand_std.mag(v,1) = mag_std;
        FRF_data.(T(n)).grand_std.gain(v,1) = gain_std;
        FRF_data.(T(n)).grand_std.phase(v,1) = phase_std;
        FRF_data.(T(n)).grand_std.time_diff(v,1) = time_diff_std;
        FRF_data.(T(n)).grand_std.coherence(:,v) = cohr_std;
    	FRF_data.(T(n)).grand_std.IO_coherence(v,1) = IOcohr_std;
        FRF_data.(T(n)).grand_std.complex(v,1) = std(FRF_data.(T(n)).fly.complex(v,:));
        FRF_data.(T(n)).grand_std.error(v,1) = std(FRF_data.(T(n)).fly.error(v,:));
    end
    
    % Fit time constant
    for f = 1:N.fly
        X = FRF_data.IOFv(2:end-1);
        Y = FRF_data.(T(n)).fly.phase(2:end-1,f);
        [fitresult, gof] = time_constant_fit(X, Y, false);
        FRF_data.(T(n)).fly.time_constant(1,f) = fitresult.b;
        FRF_data.(T(n)).fly.time_constant_r2(1,f) = gof.rsquare;
    end
    FRF_data.(T(n)).grand_mean.time_constant(1,f) = mean(FRF_data.(T(n)).fly.time_constant);
    FRF_data.(T(n)).grand_mean.time_constant_r2(1,f) = mean(FRF_data.(T(n)).fly.time_constant_r2);
    FRF_data.(T(n)).grand_std.time_constant(1,f) = std(FRF_data.(T(n)).fly.time_constant);
    FRF_data.(T(n)).grand_std.time_constant_r2(1,f) = std(FRF_data.(T(n)).fly.time_constant_r2);
end

%% FRF: one condition
cc = hsv(n_plot);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.2*n_plot 5*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(6,n_plot);
for n = 1:n_plot
    subI = n + (0:5)*n_plot;
    ax(1,n) = subplot(6,n_plot,subI(1)); hold on ; title(T(n), 'interpreter', 'none')
        plot(FRF_data.IOFv, FRF_data.ref_mag, '*-', 'Color', 'k', 'LineWidth', 0.5)
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.mag, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(1,n),h.line(1,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.mag,...
                  FRF_data.(T(n)).grand_std.mag, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(2,n) = subplot(6,n_plot,subI(2)); hold on
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.gain, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(2,n),h.line(2,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.gain,...
                  FRF_data.(T(n)).grand_std.gain, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(3,n) = subplot(6,n_plot,subI(3)); hold on
        plot([0 20], [0 0], '--k')
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.phase, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(3,n),h.line(3,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.phase,...
                  FRF_data.(T(n)).grand_std.phase, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);

    ax(4,n) = subplot(6,n_plot,subI(4)); hold on
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.time_diff, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.time_diff,...
                  FRF_data.(T(n)).grand_std.time_diff, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
    ax(5,n) = subplot(6,n_plot,subI(5)); hold on
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.error, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.error,...
                  FRF_data.(T(n)).grand_std.error, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
              
    ax(6,n) = subplot(6,n_plot,subI(6)); hold on
        plot(FRF_data.IOFv, FRF_data.(T(n)).fly.IO_coherence, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        [h.patch(6,n),h.line(6,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.IO_coherence,...
                  FRF_data.(T(n)).grand_std.IO_coherence, FRF_data.IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    
        %plot(FRF_data.Fv, FRF_data.(T(n)).fly.coherence, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
        %[h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(T(n)).grand_mean.coherence(2:end),...
                  %FRF_data.(T(n)).grand_std.coherence(2:end), FRF_data.Fv(2:end), 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
end

set(h.line(1:6,:), 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15')
set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.5 25],...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
linkaxes(ax(1,:),'y')
% linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
linkaxes(ax(4,:), 'y')
linkaxes(ax(5,:), 'y')
linkaxes(ax(6,:), 'y')

XLabelHC = get(ax(end,:), 'XLabel');
if n_plot == 1
    set([XLabelHC], 'String', 'Frequency (Hz)')
else
    set([XLabelHC{:}], 'String', 'Frequency (Hz)')
end

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Time difference (ms)')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Error')
YLabelHC = get(ax(6,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(1,1:end),'YLim',[0 3.2])
set(ax(2,1:end),'YLim',[0 1])
% set(ax(3,1:end),'YLim',[-500 200])
set(ax(3,1:end),'YLim',[-550 200])
% set(ax(4,1:end),'YLim',400*[-1 1])
set(ax(5,1:end),'YLim',[0 1.5])
set(ax(6,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])

set(ax,'XScale','log')
% set(ax,'XScale','linear')
align_Ylabels(fig)

%% FRF: compare
clc

trf_names = ["ref2body", "ref2gaze", "ref2head",];
n_plot = length(trf_names);
% cc = hsv(n_plot);
cc = [0.9 0 0 ; 0.5 0.3 1 ; 0 0.4 1 ; 0 0.8 0.2];


fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*1 5*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(5,1);
ax(1,1) = subplot(5,1,1); hold on
    for n = 1:n_plot
        [h.patch(1,n),h.line(1,n)] = PlotPatch(FRF_data.(trf_names(n)).grand_mean.mag,...
                  FRF_data.(trf_names(n)).grand_std.mag, FRF_data.IOFv, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end

ax(2,1) = subplot(5,1,2); hold on
    for n = 1:n_plot
        [h.patch(2,n),h.line(2,n)] = PlotPatch(FRF_data.(trf_names(n)).grand_mean.gain,...
                  FRF_data.(trf_names(n)).grand_std.gain, FRF_data.IOFv, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    
ax(3,1) = subplot(5,1,3); hold on
    yline(0, '--k', 'LineWidth', 1)
    for n = 1:n_plot
        [h.patch(3,n),h.line(3,n)] = PlotPatch(FRF_data.(trf_names(n)).grand_mean.phase,...
                  FRF_data.(trf_names(n)).grand_std.phase, FRF_data.IOFv, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    
ax(4,1) = subplot(5,1,4); hold on
%     ylabel('Error')
    yline(1, '--k', 'LineWidth', 1)
    for n = 1:n_plot
        [h.patch(4,n),h.line(4,n)] = PlotPatch(FRF_data.(trf_names(n)).grand_mean.error,...
                  FRF_data.(trf_names(n)).grand_std.error, FRF_data.IOFv, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end

ax(5,1) = subplot(5,1,5); hold on
    for n = 1:n_plot
        [h.patch(5,n),h.line(5,n)] = PlotPatch(FRF_data.(trf_names(n)).grand_mean.IO_coherence,...
                  FRF_data.(trf_names(n)).grand_std.IO_coherence, FRF_data.IOFv, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    

leg_label = trf_names;
leg = legend(squeeze(h.line(1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
leg.Position  = [0.35 0.93 0.51 0.055];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end

delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 11, 'LineWidth', 1.5)
set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'XLim', [0.1 10],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])
set(ax,'XScale','log')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Magnitude (°/s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', 'Error')
YLabelHC = get(ax(5,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(1,1:end),'YLim',[0 65])
set(ax(2,1:end),'YLim',[0 1])
set(ax(3,1:end),'YLim',[-200 200])
set(ax(4,1:end),'YLim',[0 1.5])
set(ax(5,1:end),'YLim',[0 1])
set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])

set(ax(1:end-1,:), 'XColor', 'none')

% align_Ylabels(fig)

%% Complex
clc

T = ["ref2body", "ref2head","ref2gaze"];
pI = [1 2 3];
yL = [1 1 1];

% T = ["ref2head", "ref2wing"];
% pI = [1 2];
% yL = [1 1];

n_plot = length(pI);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 3])
mrksz = 6;

movegui(fig, 'center')
clear ax h hh
n_freq = length(FRF_data.IOFv);
ax = gobjects(1, n_plot);
h = gobjects(n_freq, n_plot);
for n = 1:n_plot
        R = real(FRF_data.(T(n)).fly.complex);
        I = imag(FRF_data.(T(n)).fly.complex);
        [n_freq, n_fly] = size(FRF_data.(T(n)).fly.complex);
        G = repmat(FRF_data.IOFv, [n_fly, 1]);
        cc = hsv(n_freq);
        
        ax(n) = subplot(1,n_plot,n);
        [ax_cm,h_cm] = ComplexAxes(0:0.2:yL(n));
        set(ax_cm, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
        plot(1, 0, 'g.',  'MarkerSize', 25, 'MarkerFaceColor', 'none')
        
        h(:,n) = gscatter(R(:), I(:), G, cc, '.', 7, false);
        mu = FRF_data.(T(n)).grand_mean.complex;
        plot(real(mu), imag(mu), '-k.', 'MarkerSize', 12, 'MarkerFaceColor', 'none', 'LineWidth', 1);
            
        for f = 1:n_freq          
            x = [0 real(mu(f))];
            y = [0 real(mu(f))];
            theta = rad2deg(atan2(x(2), y(2)));
            gain_mu = FRF_data.(T(n)).grand_mean.gain(f);
            gain_std = FRF_data.(T(n)).grand_std.gain(f);
            phase_std = FRF_data.(T(n)).grand_std.phase(f);
            b = 2*pi*gain_mu*(phase_std/ 360);
            n_std = 1;
            elps = Ellipse([real(mu(f)) imag(mu(f))], 2*n_std*gain_std, 2*n_std*b, ...
                'ellipse', theta, 0, cc(f,:), 0.4);

            elps = draw(elps);
            set(elps.h.patch, 'EdgeColor', 'k')
            set(elps.h.centroid, 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
            hh.patch(v,n,f) = elps.h.patch;
            hh.med(v,n,f) = elps.h.centroid;
        end
        if n==1
            leg = legend(squeeze(h(:,n)), string(FRF_data.IOFv), 'Box', 'off', 'Location', 'east');
            leg.Position(1:2) = leg.Position(1:2) + [-0.2 0.04];
            leg.Title.String = 'Frequency (hz)';
        end
        title(T(n))
%         uistack(squeeze(hh.patch(v,n,:)), 'top')
%         uistack(squeeze(h(v,n,:)), 'top')
end

%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'Color', 'none', 'LineWidth', 1.5)

%% Time constant
time_const = [];
r2 = [];
time_const.body_head_free = 1000*FRF_data.(T(1)).fly.time_constant';
time_const.gaze_head_free = 1000*FRF_data.(T(3)).fly.time_constant';

r2.body_head_free = FRF_data.(T(1)).fly.time_constant_r2';
r2.gaze_head_free = FRF_data.(T(3)).fly.time_constant_r2';

time_const_all = [time_const.body_head_free ; time_const.gaze_head_free ];
r2_all = [r2.body_head_free ; r2.gaze_head_free ];
G = [1*ones(size(time_const.body_head_free)) ; 2*ones(size(time_const.gaze_head_free))];

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 4 2*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,1);
ax(1) = subplot(2,1,1); hold on
    b = boxchart(time_const_all,'GroupByColor', G, 'MarkerStyle', '.');

ax(2) = subplot(2,1,2); hold on
    b = boxchart(r2_all,'GroupByColor', G, 'MarkerStyle', '.');
fnames = fieldnames(time_const);
leg = legend(fnames, 'Box', 'off', 'interpreter', 'none', 'Orientation', 'horizontal');
leg.Position = [0.23 0.95 0.58 0.05];

linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end
set(ax(1,:), 'YLim', [-80 0])
set(ax(2,:), 'YLim', [0 1])

set(ax, 'Color', 'none', 'LineWidth', 1.2, 'FontSize', 10, 'Box', 'off', 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Time constant (ms)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'R^{2}')

%% Save FRF data
filedata = textscan(FILE, '%s', 'delimiter', '_');
dataset_name = [];
for n = 1:4
    dataset_name = [dataset_name '_' char(filedata{1}(n))];
end
fname = ['FRF' dataset_name];
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'FRF_data', 'FUNC', 'U', 'N');

end