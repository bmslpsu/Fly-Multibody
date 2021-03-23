function [] = fitTF(trfn)
%% fitTF:
%   trfn: name of transfrom (ref2body, etc.)
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FRF_data','FUNC','U','N');

%% Fit TF to mean fly frequency response functions
clearvars -except PATH FILE FRF_data FUNC U N
clc

trfn = ["ref2body"];
% trfn = ["ref2head"];
% trfn = ["ref2gaze"];
n_cond = N{1,3};
n_fly = N.fly;

if strcmp(trfn, "ref2head")
    init_sys = idtf(NaN(1,3), [1,NaN(1,2)]); % head
    phase_lim = [200 100];
    m = 0.001;
    wn = 100;
    zeta = 1;
    P = tf([1/m 0], [1 2*zeta*wn wn^(2)]);
elseif strcmp(trfn, "ref2body")
    init_sys = idtf(NaN(1,3), [1,NaN(1,2)]); % body
  	phase_lim = [0 0];
    m = 0.001;
    tau = 0.01;
    %P = tf(1/m, [tau 1]);
    P = tf(1/m, [tau 1 0]);
elseif strcmp(trfn, "ref2gaze")
    %init_sys = idtf(NaN(1,2), [1,NaN(1,3)]); % gaze
end
% init_sys = idtf(NaN(1,2), [1,NaN(1,2)], 'IODelay', NaN);
% init_sys.Structure(1).IODelay.Free = true;
% init_sys.Structure(1).IODelay.Maximum = 0.1;

reg.Lambda = 0;
reg.R = 1;
reg.Nominal = 'model';
opt  = tfestOptions('EnforceStability', true, 'Regularization', reg); 

Fv_fit = (0.1:0.1:30)';
TF_data = [];
for v = 1:n_cond
    for f = 1:n_fly
        % Fit TF
        response = FRF_data.(trfn(1)).fly(v).complex(:,f);
        sys = frd(response, FRF_data.IOFv{v}, 'FrequencyUnit', 'hz');
        
        %[b,a] = invfreqs(response, FRF_data.IOFv{v}, 2, 2);
        
        tffit = tfest(sys, init_sys, opt);
        %tffit = tf(b, a);
        if length(tffit.Denominator) == length(init_sys.Denominator)
            TF_data.fly{v}.(trfn(1))(f).sys = tffit;
            temp_fit = tffit.Report.Fit.FitPercent;
        else
            TF_data.fly{v}.(trfn(1))(f).sys = init_sys;
            temp_fit = nan;
        end
        
        % Bode gain & phase
        fly_tf = TF_data.fly{v}.(trfn(1))(f).sys;
        [gain,phase,~] = bode(fly_tf, 2*pi*Fv_fit);
        
        phase = squeeze(phase);
        phase_shift = linspace(phase_lim(1), phase_lim(2), length(phase))';
        phase(phase > phase_shift) = phase(phase > phase_shift) - 360;
        
        TF_data.fly{v}.(trfn(1))(f).gain = squeeze(gain);
        TF_data.fly{v}.(trfn(1))(f).phase = phase;
        TF_data.fly{v}.(trfn(1))(f).delay = fly_tf.IODelay;
        TF_data.fly{v}.(trfn(1))(f).num = fly_tf.Numerator';
        TF_data.fly{v}.(trfn(1))(f).den = fly_tf.Denominator';
        TF_data.fly{v}.(trfn(1))(f).fit = temp_fit;
        
        % Controller
      	fly_controller = minreal(fly_tf / ((1 - fly_tf)*P));
        [gain,phase,~] = bode(fly_controller, 2*pi*Fv_fit);
        phase = squeeze(phase);
        phase_shift = linspace(phase_lim(1), phase_lim(2), length(phase))';
        phase(phase > phase_shift) = phase(phase > phase_shift) - 360;
        TF_data.fly{v}.(trfn(1))(f).controller = fly_controller;        
        TF_data.fly{v}.(trfn(1))(f).controller_gain = squeeze(gain);
        TF_data.fly{v}.(trfn(1))(f).controller_phase = squeeze(phase);
    end
    
    % Collect fly fits
    fnames = fields(TF_data.fly{v}.(trfn(1)));
    for n = 1:length(fnames)
        TF_data.fly_all.(trfn(1))(v).(fnames{n}) = ...
            cat(2, TF_data.fly{v}.(trfn(1))(:).(fnames{n}));
        if ~strcmp(fnames{n}, 'sys') && ~strcmp(fnames{n}, 'controller')
            if strcmp(fnames{n}, 'phase')
                TF_data.fly_stats.(trfn(1))(v).(fnames{n}) = ...
                    system_stats(deg2rad(TF_data.fly_all.(trfn(1))(v).(fnames{n})), 2);
            else
                TF_data.fly_stats.(trfn(1))(v).(fnames{n}) = ...
                    system_stats(TF_data.fly_all.(trfn(1))(v).(fnames{n}), 2);
            end
        end
    end
    
    % Get median fit from TF coefficents
    b = TF_data.fly_stats.(trfn(1))(v).num.median';
    a = TF_data.fly_stats.(trfn(1))(v).den.median';
    iodelay = TF_data.fly_stats.(trfn(1))(v).delay.median';
    TF_data.grand_coeff.(trfn(1))(v).sys = tf(b, a, 'IODelay', iodelay);
    
    % Bode gain & phase for median coefficents fit
    grand_coeff_tf = TF_data.grand_coeff.(trfn(1))(v).sys;
    [gain,phase,~] = bode(grand_coeff_tf, 2*pi*Fv_fit);
    
    phase = squeeze(phase);
    phase_shift = linspace(phase_lim(1), phase_lim(2), length(phase))';
    phase(phase > phase_shift) = phase(phase > phase_shift) - 360;
        
    TF_data.grand_coeff.(trfn(1))(v).gain = squeeze(gain);
    TF_data.grand_coeff.(trfn(1))(v).phase = phase;
    TF_data.grand_coeff.(trfn(1))(v).delay = grand_coeff_tf.IODelay;
    TF_data.grand_coeff.(trfn(1))(v).num = grand_coeff_tf.Numerator{1}';
    TF_data.grand_coeff.(trfn(1))(v).den = grand_coeff_tf.Denominator{1}';
    
    % Fit TF to grand mean response
    response = FRF_data.(trfn(1)).grand_mean(v).complex;
    sys = frd(response,FRF_data.IOFv{v}, 'FrequencyUnit', 'hz');
    TF_data.grand.(trfn(1))(v).sys = tfest(sys, init_sys, opt);
    
    % Bode gain & phase for median coefficents fit
    grand_tf = TF_data.grand.(trfn(1))(v).sys;
    [gain,phase,~] = bode(grand_tf, 2*pi*Fv_fit);
    
    phase = squeeze(phase);
    phase_shift = linspace(phase_lim(1), phase_lim(2), length(phase))';
    phase(phase > phase_shift) = phase(phase > phase_shift) - 360;
        
    TF_data.grand.(trfn(1))(v).gain = squeeze(gain);
    TF_data.grand.(trfn(1))(v).phase = phase;
    TF_data.grand.(trfn(1))(v).delay = grand_tf.IODelay;
    TF_data.grand.(trfn(1))(v).num = grand_tf.Numerator';
    TF_data.grand.(trfn(1))(v).den = grand_tf.Denominator';
    TF_data.grand.(trfn(1))(v).fit = grand_tf.Report.Fit.FitPercent;
    
    % Controller
    grand_controller = minreal(grand_tf / ((1 - grand_tf)*P));
    [gain,phase,~] = bode(grand_controller, 2*pi*Fv_fit);
    phase = squeeze(phase);
    phase_shift = linspace(phase_lim(1), phase_lim(2), length(phase))';
    phase(phase > phase_shift) = phase(phase > phase_shift) - 360;
    TF_data.grand.(trfn(1))(v).controller = grand_controller;        
    TF_data.grand.(trfn(1))(v).controller_gain = squeeze(gain);
    TF_data.grand.(trfn(1))(v).controller_phase = squeeze(phase);
    
    disp(['For vel = ' num2str(U.vel{1}(v)) ])
    grand_tf
end

%% Grand response
cc_fly = distinguishable_colors(N.fly);
cc_fly = [cc_fly , 0.5*ones(size(cc_fly,1),1)];
cc_data = [0 0 0];

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 2*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,n_cond);
for v = 1:n_cond
    subI = v + (0:2)*n_cond;
    ax(1,v) = subplot(2,n_cond,subI(1)); hold on
        title([num2str(U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
        %plot(HeadFree.FRF_data.IOFv{v}, HeadFree.FRF_data.(trf_names(1)).fly(v).gain, 'Color', [1 0.5 0.5 1])
        h.fly = plot(Fv_fit, TF_data.fly_all.(trfn(1))(v).gain, 'LineWidth', 0.25);
        set(h.fly, {'Color'}, num2cell(cc_fly,2))
        [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(FRF_data.(trfn(1)).grand_mean(v).gain,...
                  FRF_data.(trfn(1)).grand_std(v).gain, FRF_data.IOFv{v}, ...
                  1, 1, cc_data, 0.7*cc_data, 0.2, 1);
        h.grand(1,v) = plot(Fv_fit, TF_data.grand.(trfn(1))(v).gain, 'r', 'LineWidth', 1.5);
        
    
    ax(2,v) = subplot(2,n_cond,subI(2)); hold on
        yline(0, '--k', 'LineWidth', 0.5)
        h.fly = plot(Fv_fit, TF_data.fly_all.(trfn(1))(v).phase, 'LineWidth', 0.25);
        set(h.fly, {'Color'}, num2cell(cc_fly,2))

        [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(FRF_data.(trfn(1)).grand_mean(v).phase,...
                  FRF_data.(trfn(1)).grand_std(v).phase, FRF_data.IOFv{v}, ...
                  1, 1, cc_data, 0.7*cc_data, 0.2, 1);
        h.grand(2,v) = plot(Fv_fit, TF_data.grand.(trfn(1))(v).phase, 'r', 'LineWidth', 1.5);
end
linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end
% delete(h.patch)

set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.1 30],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

set(ax(1,1:end),'YLim',[0 1])
set(ax(2,1:end),'YLim',[-200 200])

set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])
set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')
% set(ax,'XScale','linear')
% align_Ylabels(fig)

%% Controller response
cc_fly = distinguishable_colors(N.fly);
cc_fly = [cc_fly , 0.5*ones(size(cc_fly,1),1)];
cc_data = [0 0 0];

fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2.5*n_cond 2*2])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,n_cond);
for v = 1:n_cond
    subI = v + (0:2)*n_cond;
    ax(1,v) = subplot(2,n_cond,subI(1)); hold on
        title([num2str(U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
        %plot(HeadFree.FRF_data.IOFv{v}, HeadFree.FRF_data.(trf_names(1)).fly(v).gain, 'Color', [1 0.5 0.5 1])
        h.fly = plot(Fv_fit, TF_data.fly_all.(trfn(1))(v).controller_gain, 'LineWidth', 0.25);
        set(h.fly, {'Color'}, num2cell(cc_fly,2))
%         [h.patch(1,v,1),h.line(1,v,1)] = PlotPatch(FRF_data.(trfn(1)).grand_mean(v).gain,...
%                   FRF_data.(trfn(1)).grand_std(v).gain, FRF_data.IOFv{v}, ...
%                   1, 1, cc_data, 0.7*cc_data, 0.2, 1);
        h.grand(1,v) = plot(Fv_fit, TF_data.grand.(trfn(1))(v).controller_gain, 'r', 'LineWidth', 1.5);
        
    
    ax(2,v) = subplot(2,n_cond,subI(2)); hold on
        yline(0, '--k', 'LineWidth', 0.5)
        h.fly = plot(Fv_fit, TF_data.fly_all.(trfn(1))(v).controller_phase, 'LineWidth', 0.25);
        set(h.fly, {'Color'}, num2cell(cc_fly,2))

%         [h.patch(2,v,1),h.line(2,v,1)] = PlotPatch(FRF_data.(trfn(1)).grand_mean(v).phase,...
%                   FRF_data.(trfn(1)).grand_std(v).phase, FRF_data.IOFv{v}, ...
%                   1, 1, cc_data, 0.7*cc_data, 0.2, 1);
        h.grand(2,v) = plot(Fv_fit, TF_data.grand.(trfn(1))(v).controller_phase, 'r', 'LineWidth', 1.5);
end
linkaxes(ax, 'x')
for a = 1:size(ax,1)
    linkaxes(ax(a,:), 'y')
end
% delete(h.patch)

% set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 15, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'XLim', [0.08 30],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase difference (°)')

% set(ax(1,1:end),'YLim',[0 0.005])
set(ax(2,1:end),'YLim',[-200 200])

set(ax(1:end-1,:), 'XTickLabel', [])
set(ax(:,2:end), 'YTickLabels', [])
set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')

%% Fit pertentages
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 2])
movegui(fig, 'center')
clear ax h bx
ax(1) = subplot(1,1,1); hold on ; ylabel('Fit percent')
    fit_data = cat(1,TF_data.fly_all.(trfn(1))(:).fit)';
    fit_data = fit_data(:);
    bx(1) = boxchart(fit_data);
    
% ax(2) = subplot(2,1,2); hold on ; ylabel('Delay (ms)')
%     delay_data = 1000*cat(1,TF_data.fly_all.(trfn(1))(:).delay)';
%     delay_data = delay_data(:);
%     bx(2) = boxchart(delay_data);

linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 10, 'Box', 'off', 'XColor', 'w')
set(ax(1), 'YLim', [0 100])
% set(ax(2), 'YLim', [0 100])

%% Find controller
clc
bopt = bodeoptions;
% set(bopt, 'FreqUnits', 'Hz', 'MagUnits', 'abs', 'XLim', [0.2 20], ...
%     'PhaseWrappingBranch', -20, 'PhaseWrapping', 'on')
set(bopt, 'FreqUnits', 'Hz', 'MagUnits', 'abs', 'XLim', [0.2 20], ...
    'PhaseWrappingBranch', -200, 'PhaseWrapping', 'on')
% plant_order = [2 0];
% dc_num = 1;
% control_type = 'PID';
% clsys = TF_data.grand.ref2head(1).sys;
% [K] = get_PID(clsys, plant_order, dc_num, control_type);
% syms m tau s
m = 0.001;
tau = 0.01;
wn = 100;
zeta = 1;
P = tf(1/m, [tau 1]);
% P = tf(1/m, [1 2*zeta*wn wn^(2)]);
G = (TF_data.grand.(trfn(1))(3).sys);
% P = symtf(sym(1/m), sym([tau 1]));
% G = vpa(collect(tf2sym(tf(TF_data.grand.(trfn(1))(v).sys))),3);
C = G / ((1 - G)*P);
C = minreal(C);
% C = vpa(collect(simplify(expand(C)), s), 3);
% C = sym2tf(C, s)

close all
bodeplot(C, bopt)

get_controller(G, P)

%% Stats
% [p,tb,stats] = anova1(time_const_keep, G_keep);
% % [p,tb,stats] = kruskalwallis(time_const_keep, G_keep);
% [c,m] = multcompare(stats, 'alpha', 0.001);
 
%% Save TF fit data
filedata = textscan(FILE, '%s', 'delimiter', '._');
filename = [];
for n = 2:5
    filename = [filename '_' char(filedata{1}(n))];
end
fname = ['TFfit_' char(trfn(1)) filename]; 

root = 'E:\DATA\Magno_Data\Multibody';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'TF_data','FUNC','U','N');
end