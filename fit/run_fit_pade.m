function [] = run_fit_pade()
%% run_fit_pade:
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
showplot = true;

%% HeadFree: err2body
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'err2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 0, 0.02, opt);
tffit{end+1} = tfest(data, 1, 1, 0.02, opt);

[model, h] = fit_complex(Cn, IOFv, [1 0], [1 1], 0.02, 'GP', 'lsqcurvefit', false);
tffit{end+1} = model.G;
% [model, h] = fit_complex(Cn, IOFv, 1, [1 1], 0.02, 'GP', 'lsqcurvefit', true);
% [model, h] = fit_tf(Cn, IOFv, [1], [1 1], 0.02, 'GP', 'lsqnonlin', true);

% delay = 0:0.001:0.05;
% [tffit{end+1}, sys_list, fitpercent, delay_sort] = tfest_delay(data, 1, 0, opt, delay);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    suptitle(trf)
    h.fig.Position(3:4) = [2.5 6];
    set(h.data, 'MarkerSize', 5)
    set(h.sys_freq, 'MarkerSize', 10)
    set(h.ax, 'LineWidth', 0.5)
    % set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
    set(h.ax(1), 'YLim', [-7 1], 'XLim', [-2 6])
    set(h.ax(2), 'YLim', [0 5])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: err2head
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'err2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
% delay = 0:0.001:0.05;
% [tffit{end+1}, sys_list, fitpercent, delay_sort] = tfest_delay(data, 1, 1, opt, delay);

tffit{end+1} = tfest(data, 1, 1, 0.02, opt);
tffit{end+1} = tfest(data, 1, 1, 0.02, opt);
tffit{end+1} = tfest(data, 1, 0, 0.02, opt);
tffit{end+1} = tfest(data, 2, 1, 0.02, opt);
% tffit{1} = tf([tffit{1}.Numerator(1) 0], tffit{1}.Denominator, 'IODelay', 0.02);
tffit{1}.Numerator = [tffit{1}.Numerator(1) 0];
% tffit{3}.Numerator = [tffit{3}.Numerator(1) 0];
% tffit{end+1} = tfest(data, 2, 1, 0.02, opt);
% tffit{end+1} = tfest(data, 2, 2, [], opt);

% [temp, h] = fitBode(Cn, IOFv, [1 1], [1 1], 0.02, 'lsqnonlin', false);
% tffit{end+1} = temp.tf.G;

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    suptitle(trf)
    h.fig.Position(3:4) = [3 7];
    set(h.data, 'MarkerSize', 5)
    set(h.sys_freq, 'MarkerSize', 10)
    set(h.ax, 'LineWidth', 0.5)
    % set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: ref2body
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 0, 0.02, opt);
% tffit{end+1} = tfest(data, 1, 0, 0.02, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: ref2head
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 1, 0.02, opt);
tffit{1}.Numerator = [tffit{1}.Numerator(1) 0];
% tffit{end+1} = tfest(data, 1, 2, [], opt);
% tffit{end+1} = tfest(data, 2, 1, [], opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: ref2gaze
close all
clc
clear Cn sys
clss = 'HeadFree';
trf = 'ref2gaze';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 2, 2, 0.02, opt);
% tffit{end+1} = tfest(sys, 1, 2, NaN, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFixed: err2body
close all
clc
clear Cn sys
clss = 'HeadFixed';
trf = 'err2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];

% delay = 0:0.001:0.04;
% [best_fit, sys_list, fitpercent, delay_sort] = tfest_delay(data, 1, 0, opt, delay);

tffit{end+1} = tfest(data, 1, 0, 0.023, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    % set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 3])
    set(h.ax(3), 'YLim', [-250 150])
end

%% BodyFixed: err2head
close all
clc
clear Cn sys
clss = 'BodyFixed';
trf = 'err2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];

% delay = 0:0.001:0.05;
% [best_fit, sys_list, fitpercent, delay_sort] = tfest_delay(data, 1, 1, opt, delay);

tffit{end+1} = tfest(data, 1, 1, 0.029, opt);
% tffit{end+1} = tfest(data, 2, 1, 0.032, opt);
tffit{1}.Numerator = [tffit{1}.Numerator(1) 0];
% tffit{1}.Denominator = [1 8];
% tffit{end+1} = tfest(sys, 1, 1, NaN, opt);
% pp = pole(tffit{1});
% 
% tffit{end+1} = tf([tffit{1}.Numerator(1) 0], [1 6]);
% tffit{end}.IODelay = tffit{1}.IODelay;


% delay = 0:0.001:0.05;
% [tffit{end+1}, sys_list, fitpercent, delay_sort] = tfest_delay(data, 2, 1, opt, delay);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    % set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 2])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFixed: ref2body
close all
clc
clear Cn sys
clss = 'HeadFixed';
trf = 'ref2body';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 0, 0.02, opt);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    % set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% BodyFixed: ref2head
close all
clc
clear Cn sys
clss = 'BodyFixed';
trf = 'ref2head';
Cn = {ALL.(clss).FRF_data.(trf).fly(vI).gain , ALL.(clss).FRF_data.(trf).fly(vI).phase};
data = frd(ALL.(clss).FRF_data.(trf).grand_mean(vI).complex, IOFv, 'FrequencyUnit', 'hz');

tffit = [];
tffit{end+1} = tfest(data, 1, 1, 0.02, opt);
tffit{1}.Numerator = [tffit{1}.Numerator(1) 0];
% delay = 0:0.001:0.05;
% [tffit{end+1}, sys_list, fitpercent, delay_sort] = tfest_delay(data, 1, 1, opt, delay);

[MODEL.(clss).fit.(trf), MODEL.(clss).data.(trf), h] = plotFit(Cn, IOFv, tffit, frange, showplot);

if showplot
    % set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% Create closed-loop transforms from open-loop fits
pN = 1;
clss = 'HeadFree';
MODEL.morph.body.M = 0.84185066; % [mg] body mass
MODEL.morph.head.M = 0.08957730; % [mg] head mass
MODEL.morph.body.L = 2.22; % [mm] body length
MODEL.morph.head.L = 0.42; % [mm] head length
MODEL.morph.body.R = 0.8; % [mm] body radius
MODEL.morph.head.R = 0.65; % [mm] head radius
MODEL.morph.body.Jzz = 0.56650720*(1e-3)^(3); % [kg*m^2] body inertia
MODEL.morph.head.Jzz = 0.00619410*(1e-3)^(3); % [kg*m^2] head inertia

MODEL.morph.body.Jzz_fromMass = MODEL.morph.body.M*( (1/12)*(MODEL.morph.body.L)^(2) + ...
                                    (MODEL.morph.body.L/3)^(2) + (1/4)*MODEL.morph.body.R^(2)) * (1e-3)^(3); % [kg*m^2] body inertia
                                
MODEL.morph.head.Jzz_fromMass = MODEL.morph.head.M*( (1/12)*(MODEL.morph.head.L)^(2) + ...
                                    0.2*(MODEL.morph.body.L)^(2) + (1/4)*MODEL.morph.head.R^(2)) * (1e-3)^(3); % [kg*m^2] head inertia

MODEL.morph.J_ratio = MODEL.morph.head.Jzz / MODEL.morph.body.Jzz; % head/body inertia ratio

MODEL.(clss).G.body = MODEL.(clss).fit.err2body(1).models;
MODEL.(clss).G.body_bad = MODEL.(clss).fit.err2body(2).models;
MODEL.(clss).G.head = MODEL.(clss).fit.err2head(1).models;
MODEL.(clss).G.head_bad = MODEL.(clss).fit.err2head(2).models;
MODEL.(clss).G.head_2nd = MODEL.(clss).fit.err2head(3).models;
MODEL.(clss).G.gaze = MODEL.(clss).G.body + MODEL.(clss).G.head;

MODEL.(clss).G.pade.body = pade(MODEL.(clss).G.body, pN);
MODEL.(clss).G.pade.head = pade(MODEL.(clss).G.head, pN);

[MODEL.(clss).C.body, MODEL.(clss).P.body] = get_controller(MODEL.(clss).G.body, nan);
[MODEL.(clss).C.head, MODEL.(clss).P.head] = get_controller(MODEL.(clss).G.head, nan);

[MODEL.(clss).C_norm.head, MODEL.(clss).P_norm.head] = get_controller(MODEL.(clss).G.head, MODEL.morph.head.Jzz);
[MODEL.(clss).C_norm.body, MODEL.(clss).P_norm.body] = get_controller(MODEL.(clss).G.body, MODEL.morph.body.Jzz);

MODEL.(clss).C_norm.body = tf(MODEL.(clss).G.body.numerator / MODEL.(clss).G.body.denominator(end), 1);
MODEL.(clss).C_norm.head = tf(MODEL.(clss).G.head.numerator / MODEL.(clss).G.head.denominator(end), 1);

MODEL.(clss).H.body = minreal( MODEL.(clss).G.body / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );
MODEL.(clss).H.head = minreal( MODEL.(clss).G.head / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );
MODEL.(clss).H.gaze = minreal( (MODEL.(clss).G.body + MODEL.(clss).G.head) / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );

MODEL.(clss).H.direct.body = MODEL.(clss).fit.ref2body(1).models;
MODEL.(clss).H.direct.head = MODEL.(clss).fit.ref2head(1).models;
MODEL.(clss).H.direct.gaze = MODEL.(clss).fit.ref2gaze(1).models;

MODEL.(clss).H.pade.body = minreal( MODEL.(clss).G.pade.body / (1 + MODEL.(clss).G.pade.body + MODEL.(clss).G.pade.head) );
MODEL.(clss).H.pade.head = minreal( MODEL.(clss).G.pade.head / (1 + MODEL.(clss).G.pade.body + MODEL.(clss).G.pade.head) );
MODEL.(clss).H.pade.gaze = minreal( (MODEL.(clss).G.pade.body + MODEL.(clss).G.pade.head) ...
                                    / (1 + MODEL.(clss).G.pade.body + MODEL.(clss).G.pade.head) );

MODEL.(clss).H.body_no_head = minreal( MODEL.(clss).G.body / (1 + MODEL.(clss).G.body));
MODEL.(clss).H.head_no_body = minreal( MODEL.(clss).G.head / (1 + MODEL.(clss).G.head));

MODEL.(clss).H.pade.body_no_head = minreal( MODEL.(clss).G.pade.body / (1 + MODEL.(clss).G.pade.body));
MODEL.(clss).H.pade.head_no_body = minreal( MODEL.(clss).G.pade.head / (1 + MODEL.(clss).G.pade.head));

MODEL.(clss).W.body = minreal( MODEL.(clss).C.body / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );
MODEL.(clss).W.head = ( MODEL.(clss).C.head / (1 + MODEL.(clss).G.body + MODEL.(clss).G.head) );

clss = 'HeadFixed';
MODEL.(clss).G.body = MODEL.(clss).fit.err2body(1).models;
MODEL.(clss).H.body = minreal( MODEL.(clss).G.body / (1 + MODEL.(clss).G.body));

MODEL.(clss).G.pade.body = pade(MODEL.(clss).G.body, pN);
MODEL.(clss).H.pade.body = minreal( MODEL.(clss).G.pade.body / (1 + MODEL.(clss).G.pade.body));

[MODEL.(clss).C.body, MODEL.(clss).P.body] = get_controller(MODEL.(clss).G.body, MODEL.morph.body.Jzz);

[MODEL.(clss).C.body, MODEL.(clss).P.body] = get_controller(MODEL.(clss).G.body, nan);

[MODEL.(clss).C_norm.body, MODEL.(clss).P_norm.body] = get_controller(MODEL.(clss).G.body, MODEL.morph.body.Jzz);

MODEL.(clss).H.direct.body = MODEL.(clss).fit.ref2body(1).models;

clss = 'BodyFixed';
MODEL.(clss).G.head = MODEL.(clss).fit.err2head(1).models;
MODEL.(clss).H.head = minreal( MODEL.(clss).G.head / (1 + MODEL.(clss).G.head));

MODEL.(clss).G.pade.head = pade(MODEL.(clss).G.head, pN);
MODEL.(clss).H.pade.head = minreal( MODEL.(clss).G.pade.head / (1 + MODEL.(clss).G.pade.head));

[MODEL.(clss).C.head, MODEL.(clss).P.head] = get_controller(MODEL.(clss).G.head, MODEL.morph.head.Jzz);

[MODEL.(clss).C.head, MODEL.(clss).P.head] = get_controller(MODEL.(clss).G.head, nan);

[MODEL.(clss).C_norm.head, MODEL.(clss).P_norm.head] = get_controller(MODEL.(clss).G.head, MODEL.morph.head.Jzz);

MODEL.(clss).H.direct.head = MODEL.(clss).fit.ref2head(1).models;

%% HeadFree: ref2body compare
close all
clc
trf = 'ref2body';
% models = {MODEL.HeadFree.fit.(trf)(1).models, MODEL.HeadFree.H.body};
% models = {MODEL.HeadFree.H.body, MODEL.HeadFree.H.body_no_head};
models = {MODEL.HeadFree.H.body, MODEL.HeadFree.H.body_no_head, MODEL.HeadFixed.H.body};
[sys,~,h] = plotFit(MODEL.HeadFixed.data.(trf).input, IOFv, models, frange, 1);

fitpercent = [sys.fitpercent];
disp('Fit percent:')
disp([fitpercent.combined])

if showplot
    suptitle(trf)
    h.fig.Position(3:4) = [2 5];
    set(h.data, 'MarkerSize', 5, 'Marker', '.')
    set(h.sys_freq, 'MarkerSize', 7)
    set(h.ax, 'LineWidth', 0.5)
    % set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: ref2head compare
close all
clc
trf = 'ref2head';
% models = {MODEL.HeadFree.fit.(trf)(1).models, MODEL.HeadFree.H.head};
% models = {MODEL.HeadFree.H.head, MODEL.HeadFree.H.head_no_body};
% models = {MODEL.HeadFree.H.head, MODEL.HeadFree.H.head_no_body, MODEL.BodyFixed.H.head};

models = {MODEL.HeadFree.H.pade.head, MODEL.HeadFree.H.pade.head_no_body, MODEL.BodyFixed.H.head};
[sys,~,h] = plotFit(MODEL.BodyFixed.data.(trf).input, IOFv, models, frange, true);

fitpercent = [sys.fitpercent];
disp('Fit percent:')
disp([fitpercent.combined])

if showplot
    suptitle(trf)
    h.fig.Position(3:4) = [2 5];
    set(h.data, 'MarkerSize', 5, 'Marker', 'none')
    set(h.sys_freq, 'MarkerSize', 7)
    set(h.ax, 'LineWidth', 0.5)
    % set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% HeadFree: ref2gaze compare
trf = 'ref2gaze';
% models = {MODEL.HeadFree.fit.(trf)(1).models, MODEL.(clss).H.gaze};
models = {MODEL.HeadFree.H.gaze};
[~,~,h] = plotFit(MODEL.HeadFree.data.(trf).input, IOFv, models, frange, showplot);

if showplot
    suptitle(trf)
    h.fig.Position(3:4) = [2.5 4.5];
    set(h.data, 'MarkerSize', 5)
    set(h.sys_freq, 'MarkerSize', 7)
    set(h.ax, 'LineWidth', 0.5)
    % set(h.ax, 'XGrid', 'off', 'YGrid', 'off')
    set(h.ax(1), 'YLim', [-1 1], 'XLim', [-1 1])
    set(h.ax(2), 'YLim', [0 1])
    set(h.ax(3), 'YLim', [-250 150])
end

%% Save TF fit data
fname = ['MODEL_v3_' FILE]; 
root = 'E:\DATA\Magno_Data\Multibody';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'MODEL','ALL');
end