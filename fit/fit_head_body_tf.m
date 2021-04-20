function [] = fit_head_body_tf()
%% fit_head_body:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL');

%% Set initialization parameters
clearvars -except PATH FILE ALL
clc
vI = 2;
IOFv = ALL.HeadFree.FRF_data.IOFv{vI};
frange = 0:0.02:20;

opt = tfestOptions('EnforceStability', true, 'InitializeMethod', 'all');

%% HeadFree: err2body
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'err2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

% sys0 = idtf([nan],[1 nan]);
% sys0.Structure.IODelay.Value = 0.1; % initial guess
% sys0.Structure.IODelay.Maximum = 0.15; % maximum allowable value for delay 
% sys0.Structure.IODelay.Free = true; % treat delay as estimatable quantity
% tffit = [];
% tffit{end+1} = tfest(sys, sys0);

tffit = [];
tffit{end+1} = tfest(data, 1, 0, 0.02, opt);

delay = 0.002:0.002:0.08;
[best_sys, sys_list, fitpercent] = tfest_delay(data, 1, 0, opt, delay);

% tffit{end+1} = tfest(sys, 3, 0, [], opt);
% tffit{end+1} = tfest(sys, 2, 0, [], opt);
% [model, h] = fitBode(Cn, IOFv, [1], [1 1], 0.0224, 'fmincon', false);
% [model, h] = fitBode(Cn, IOFv, [1 1 1], [1 1 1], [], 'fmincon', false);
% tffit{end+1} = model.tf.G;

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

set(h.ax(1), 'YLim', [-7 1], 'XLim', [-2 6])
suptitle(trf)
h.fig.Position(3:4) = [2.5 6];
set(h.data, 'MarkerSize', 5)
set(h.sys_freq, 'MarkerSize', 7)
set(h.ax, 'LineWidth', 0.5)
% set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
set(h.ax(2), 'YLim', [0 5])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: err2head
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'err2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 2, [], opt);
% tffit{end+1} = tfest(sys, 2, 1, NaN, opt);
% [model, h] = fitBode(Cn, IOFv, [1 1], [1 1 1], 0.0224, 'fmincon', false);

% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);


%%% z = zero(tffit{1});
% p = pole(tffit{1});
% tffit{end+1} = tf([1 z], [1 max(p)], 'IODelay', tffit{1}.IODelay);
% 
% tffit{end+1} = zpk(z,p,1);
% % tffit{end+1} = tf([1 -z(1)], [1 -max(p)], 'IODelay', tffit{1}.IODelay);
% 
% temp = tffit{1};
% z = zero(temp);
% p = pole(temp);
% den_test = [1 -max(p)];
% sys1 = tf(temp.num, temp.den);
% % sys2 = tf([1 0], temp.den);
% % sys3 = tf(1, temp.den);
% sys2 = tf(temp.num, den_test);
% sys3 = tf(1, [1 -min(p)]);
% % sys2 = tf([1], den_test);
% 
% h = bodeplot(sys1, sys2, sys3);
% % p = getoptions(h); 
% % p.PhaseMatching = 'on'; 
% % p.PhaseMatchingFreq = 1; 
% % p.PhaseMatchingValue = 90;
% % setoptions(h,p)
% setoptions(h,'FreqUnits','Hz')
% xlim([0.1 25])
% % tffit{end+1} = sys1;

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

suptitle(trf)
h.fig.Position(3:4) = [2.5 4.5];
set(h.data, 'MarkerSize', 5)
set(h.sys_freq, 'MarkerSize', 7)
set(h.ax, 'LineWidth', 0.5)
% set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: ref2body
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 1, NaN, opt);
% tffit{end+1} = tfest(sys, 2, 0, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: ref2head
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 2, [], opt);
tffit{end+1} = tfest(data, 1, 2, [], opt);
tffit{end+1} = tfest(data, 2, 1, [], opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: ref2gaze
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2gaze';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 2, [], opt);
% tffit{end+1} = tfest(sys, 1, 2, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% Create closed-loop transforms from open-loop fits
MODEL.(clss).G.body = MODEL.(clss).fit.err2body(1).models;
MODEL.(clss).G.head = MODEL.(clss).fit.err2head(1).models;

MODEL.(clss).H.body = minreal( MODEL.(clss).G.body / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );
MODEL.(clss).H.head = minreal( MODEL.(clss).G.head / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );
MODEL.(clss).H.gaze = minreal( (MODEL.(clss).G.body + MODEL.(clss).G.head) / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );

%% HeadFree: ref2body compare
trf = 'ref2body';
% models = {MODEL.(clss).fit.(trf)(1).models, MODEL.(clss).H.body};
models = {MODEL.(clss).H.body};
[~,~,h] = plotFit(MODEL.(clss).data.(trf).input, IOFv, models, frange, true);

suptitle(trf)
h.fig.Position(3:4) = [2.5 4.5];
set(h.data, 'MarkerSize', 5)
set(h.sys_freq, 'MarkerSize', 7)
set(h.ax, 'LineWidth', 0.5)
% set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: ref2head compare
trf = 'ref2head';
% models = {MODEL.(clss).fit.(trf)(1).models, MODEL.(clss).H.head};
models = {MODEL.(clss).H.head};
[~,~,h] = plotFit(MODEL.(clss).data.(trf).input, IOFv, models, frange, true);

suptitle(trf)
h.fig.Position(3:4) = [2.5 4.5];
set(h.data, 'MarkerSize', 5)
set(h.sys_freq, 'MarkerSize', 7)
set(h.ax, 'LineWidth', 0.5)
% set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFree: ref2gaze compare
trf = 'ref2gaze';
% models = {MODEL.(clss).fit.(trf)(1).models, MODEL.(clss).H.gaze};
models = {MODEL.(clss).H.gaze};
[~,~,h] = plotFit(MODEL.(clss).data.(trf).input, IOFv, models, frange, true);

suptitle(trf)
h.fig.Position(3:4) = [2.5 4.5];
set(h.data, 'MarkerSize', 5)
set(h.sys_freq, 'MarkerSize', 7)
set(h.ax, 'LineWidth', 0.5)
% set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFixed: err2body
close all
clc
clear Cn sys
clss = 'HeadFixed';
trf = 'err2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 1, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 0, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

% set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 3])
set(h.ax(3), 'YLim', [-250 150])

%% BodyFixed: err2head
close all
clc
clear Cn sys
clss = 'BodyFixed';
trf = 'err2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 2, 0, opt);
tffit{end+1} = tfest(data, 2, 1, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

% set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% HeadFixed: ref2body
close all
clc
clear Cn sys
clss = 'HeadFixed';
trf = 'ref2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 1, [], opt);
% tffit{end+1} = tfest(sys, 1, 0, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

% set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%% BodyFixed: ref2head
close all
clc
clear Cn sys
clss = 'BodyFixed';
trf = 'ref2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 1, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 0, NaN, opt);
% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, true);

% set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
set(h.ax(2), 'YLim', [0 1])
set(h.ax(3), 'YLim', [-250 150])

%%
% set(h.fig, 'CurrentAxes', h.ax(1))
% plot(R, I, 'c')
% 
% set(h.fig, 'CurrentAxes', h.ax(2))
% plot(fout, gain, 'c')
% 
% set(h.fig, 'CurrentAxes', h.ax(3))
% plot(fout, phase, 'c')

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