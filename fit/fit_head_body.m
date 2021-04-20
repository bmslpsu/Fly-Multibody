function [] = fit_head_body(trfn)
%% fit_head_body:
%   trfn: name of transfrom (ref2body, etc.)
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'ALL');

%% Fit TF to mean fly frequency response functions
clearvars -except PATH FILE ALL
clc
vI = 2;
IOFv = ALL.HeadFree.FRF_data.IOFv{vI};

%% HeadFree: err2body
clc
clear Cn
% Cn = ALL.HeadFree.FRF_data.err2body.fly(vI).complex;
Cn{1} = ALL.HeadFree.FRF_data.err2body.fly(vI).gain;
Cn{2} = ALL.HeadFree.FRF_data.err2body.fly(vI).phase;
% IOFv = ALL.HeadFree.FRF_data.IOFv{vI};

num = [1];
den = [1 1];
tau = 0.021;
% tau = [];
data_method = 'GP';
fit_method = 'lsqnonlin';
[model, h] = fit_complex(Cn, IOFv, num, den, tau, data_method, fit_method, true);
model.G
model.fitpercent

h.fig.Position(3:4) = [3 8];
% set(h.ax(2), 'YLim', [0 1])
% set(h.ax(3), 'YLim', [-250 50])
set(h.ax(2:3), 'XLim', [0.2 20])

%% HeadFree: err2head
clc
Cn = ALL.HeadFree.FRF_data.err2head.fly(vI).complex;
Cn = ALL.HeadFree.FRF_data.err2head.grand_mean(vI).complex;

num = [1 1];
den = [1 1];
tau = 0.01;
data_method = 'GP';
fit_method = 'lsqnonlin';
[model, h] = fit_complex(Cn, IOFv, num, den, tau, data_method, fit_method, true);
model.G
model.fitpercent

h.fig.Position(3:4) = [3 8];
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-200 150])
set(h.ax(2:3), 'XLim', [0.2 20])
set(h.data, 'MarkerSize', 5)

%% HeadFree: ref2body
clc
Cn = ALL.HeadFree.FRF_data.ref2body.grand_mean(vI).complex;
freq_all = repmat(IOFv, [size(Cn,2) 1]);

reg.Lambda = 0;
reg.R = 1;
reg.Nominal = 'model';
opt  = tfestOptions('EnforceStability', true, 'Regularization', reg); 
init_sys = idtf(NaN(1,1), [1,NaN(1,2)]); % body

sys = frd(ALL.HeadFree.FRF_data.ref2body.grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');
tffit = tfest(sys, init_sys, opt);

num = [1];
den = [1 1 1];
tau = [];
data_method = 'GP';
fit_method = 'lsqnonlin';
[model, h] = fit_complex(Cn, IOFv, num, den, tau, data_method, fit_method, true);
model.G
model.fitpercent

h.fig.Position(3:4) = [3 8];
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-200 50])
set(h.ax(2:3), 'XLim', [0.2 20])
set(h.data, 'MarkerSize', 15)

[gain,phase,wout] = bode(tf(tffit));
[R,I] = nyquist(tf(tffit));
gain = squeeze(gain);
phase = squeeze(phase);
fout = wout ./ (2*pi);
R = squeeze(R);
I = squeeze(I);

set(h.fig, 'CurrentAxes', h.ax(1))
plot(R, I, 'c')

set(h.fig, 'CurrentAxes', h.ax(2))
plot(fout, gain, 'c')

set(h.fig, 'CurrentAxes', h.ax(3))
plot(fout, phase, 'c')

%% HeadFree: ref2head
clc
Cn = ALL.HeadFree.FRF_data.ref2head.fly(vI).complex;

num = [1 0];
den = [1 1 1];
tau = [];
data_method = 'GP';
fit_method = 'lsqnonlin';
[model, h] = fit_complex(Cn, IOFv, num, den, tau, data_method, fit_method, true);
model.G
model.fitpercent

h.fig.Position(3:4) = [3 8];
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-200 150])
set(h.ax(2:3), 'XLim', [0.2 20])
set(h.data, 'MarkerSize', 5)

%%

H = ( model.G / (1 + model.G) );

[gain,phase,wout] = bode(H);
[R,I] = nyquist(H);
gain = squeeze(gain);
phase = squeeze(phase);
fout = wout ./ (2*pi);
R = squeeze(R);
I = squeeze(I);

%% Save TF fit data
% filedata = textscan(FILE, '%s', 'delimiter', '._');
% filename = [];
% for n = 2:5
%     filename = [filename '_' char(filedata{1}(n))];
% end
% fname = ['TFfit_' char(trfn(1)) filename]; 
% 
% root = 'E:\DATA\Magno_Data\Multibody';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'TF_data','FUNC','U','N');
end