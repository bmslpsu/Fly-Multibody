function [] = fit_adaptive_model()
%% fit_adaptive_model:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

%% Set initialization parameters
clearvars -except PATH FILE ALL MODEL
clc
vI = 2;
IOFv = ALL.HeadFree.FRF_data.IOFv{vI};
wj = 2*pi*IOFv*1j;
frange = 0:0.02:20;

G_free = tf(MODEL.HeadFree.G.head);
G_free_sym = vpa( tf2sym(G_free) ,3);

% H_free = tf(MODEL.HeadFree.H.head);
% H_free_sym = vpa( tf2sym(G_free) ,3);

% syms s
% G_free_sym_freq = subs(G_free_sym, s, 2*pi*IOFv*1i);
% G_gain = abs(G_free_sym_freq);
% G_phase = angle(G_free_sym_freq);

% opt = tfestOptions('EnforceStability', true, 'InitializeMethod', 'all');
showplot = true;

%% HeadFree: err2body
close all
clc

syms s a_0 a_1 a_2
A(s) = a_0;

[N,D] = numden(G_free_sym);
nD = coeffs(D);
D = D ./ nD(2);
N = N ./ nD(2);

G_mod = exp(-s*G_free.IODelay) * (N + A) / D;
H_mod = G_mod / (1 + G_mod);
H_mod_gain_func = abs(H_mod);
H_mod_phase_func = rad2deg(angle(H_mod));

data_gain = ALL.BodyFixed.FRF_data.ref2head.grand_mean(vI).gain;
data_phase = ALL.BodyFixed.FRF_data.ref2head.grand_mean(vI).phase;
% F = [(G_mod_gain_func - data_gain).^(2) , 100*(G_mod_phase_func - data_phase).^(2)];
F = [H_mod_gain_func , H_mod_phase_func];
F = matlabFunction(F);

% fcn = @(a) sum(F(wj,a(1),a(2),a(3)));
fcn = @(a,x) F(x,a(1));

n_param = 1;
lb = -inf*ones(n_param,1)';
ub = inf*ones(n_param,1)';
x0 = ones(n_param,1)';

% opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective', 'Display', 'iter');
%[x, resnorm, residuals, exitflag, output] = lsqnonlin(objfcn, x0, lb, ub, opts);
% [Aout, resnorm, residuals, exitflag, output]  = lsqnonlin(fcn, x0, lb, ub, opts);

opts = optimoptions('lsqcurvefit', 'Display', 'iter');
X = [data_gain , data_phase];
freq_temp = repmat(2*pi*IOFv, [1, 1]);
Aout = lsqcurvefit(fcn, x0, freq_temp, X, lb, ub, opts);

% G_mod_fit = exp(-s*G_free.IODelay) * (N + Aout(1)) / D;
G_mod_fit = (N + Aout(1)) / D;
G_mod_fit = minreal(sym2tf(G_mod_fit));
G_mod_fit.IODelay = G_free.IODelay;
H_mod_fit = pade(G_mod_fit) / (1 + pade(G_mod_fit));

[gain_fit,phase_fit] = bode(H_mod_fit, 2*pi*IOFv);
gain_fit_fmin = squeeze(gain_fit);
phase_fit_fmin  = squeeze(phase_fit);

if median(phase_fit_fmin) > 150
    phase_fit_fmin = phase_fit_fmin - 360;
end

% exp(-s*G_free.IODelay)
% G_mod_freq = subs(G_mod, s, 2*pi*IOFv*1j);

% % dgain = 0.7:0.1:1.1;
% dgain = 0:0.5:10;
% Gain_fit = [];
% Phase_fit = [];
% for n = 1:length(dgain)
%     G_test = G_free;
%     G_test.Numerator{1} = [G_free.Numerator{1}(1) dgain(n)];
%     [gain_fit,phase_fit] = bode(G_test, 2*pi*IOFv);
%     Gain_fit(:,n) = squeeze(gain_fit);
%     Phase_fit(:,n) = squeeze(phase_fit);
% end

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches')
ax(1) = subplot(2,1,1) ; cla; hold on
    plot(IOFv, ALL.BodyFixed.FRF_data.ref2head.fly(vI).gain, '.-', 'Color', [0.5 0.5 0.5 0.5])
    plot(IOFv, ALL.BodyFixed.FRF_data.ref2head.grand_mean(vI).gain, '.k-', 'LineWidth', 1)
    plot(IOFv, gain_fit_fmin, 'r', 'LineWidth', 1)
    %plot(IOFv, Gain_fit,'b')
ax(2) = subplot(2,1,2) ; cla; hold on
    yline(0, '--', 'Color', [0.5 0.5 0.5]);
    plot(IOFv, ALL.BodyFixed.FRF_data.ref2head.fly(vI).phase, '.-', 'Color', [0.5 0.5 0.5 0.5])
    plot(IOFv, ALL.BodyFixed.FRF_data.ref2head.grand_mean(vI).phase, '.k-', 'LineWidth', 1)
    plot(IOFv, phase_fit_fmin, 'r', 'LineWidth', 1)
    %plot(IOFv, Phase_fit,'b')
   
linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XScale', 'log', 'XLim', [0.3 15])
set(ax(1), 'YScale', 'log')

set(ax(1)', 'YLim', [0 1.2])
set(ax(2)', 'YLim', [-200 100])


%% Save TF fit data
% fname = ['MODEL_v2_' FILE]; 
% root = 'E:\DATA\Magno_Data\Multibody';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'MODEL','ALL');
end