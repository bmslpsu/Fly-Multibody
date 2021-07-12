function [] = SOS_power_ratio()
%% SOS_power_ratio:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

%% For flies usng fit plants and controllers
clearvars -except PATH FILE ALL MODEL equations
clc
close all
syms s

% P_body = minreal(MODEL.HeadFree.P.body);
% C_body = minreal(MODEL.HeadFree.C.body);
% 
% P_head = minreal(MODEL.HeadFree.P.head);
% C_head = minreal(MODEL.HeadFree.C.head);

P_body = minreal(MODEL.HeadFree.P_norm.body);
C_body = minreal(MODEL.HeadFree.C_norm.body);

P_head = minreal(MODEL.HeadFree.P_norm.head);
C_head = minreal(MODEL.HeadFree.C_norm.head);

head_pc = C_head.Numerator{1}(1) / P_head.Numerator{1}(2);
body_tau = 1 / P_body.Denominator{1}(2);
body_static = double(subs(tf2sym(P_body),s,0));
head_static = double(subs(tf2sym(P_head),s,0));

% delay = 0.02;
% [sys.fly, data.fly, ~] = switch_system_v3(P_body, P_head, C_body, C_head, delay, delay, false);
[sys.fly, data.fly, ~] = plant_power(P_body, P_head, false);

%% Animal #1: human head-eye
close all
clc

zeta = 0.5;
wn = 2*pi*0.3;
P_head = tf(wn^(2), [1 2*zeta*wn wn^(2)]);
P_eye = tf(1, [0.13 1]);

bopt = bodeoptions; bopt.FreqUnits = 'Hz'; bopt.MagUnits = 'abs'; bopt.XLim = [0.1 20];
h = bodeplot(P_eye, P_head, bopt, 2*pi*(0:0.01:20));

[sys.human, data.human, ~] = plant_power(P_head, P_eye, true);


%% Animal comparison
close all ; clc
power_ratio = structfun(@(x) x.power.ratio_1, data, 'UniformOutput', false);
negI = structfun(@(x) sign(x), power_ratio, 'UniformOutput', false);

animal_names = fieldnames(power_ratio);
na = length(animal_names);
power_ratio_log = structfun(@(x,y) log10(abs(x)), power_ratio, 'UniformOutput', false);
for n = 1:na
    uz = power_ratio.(animal_names{n}) < 1;
    temp = power_ratio_log.(animal_names{n}) .* negI.(animal_names{n});
    temp(uz) = temp(uz) * -1;
    power_ratio_log.(animal_names{n}) = temp;
end

fig = figure (200);
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [3 4];
set(fig, 'Visible', 'on')
clear ax h

cc = hsv(na);
ax = subplot(1,1,1); hold on ; cla
for n = 1:na
    h.power(n) = plot(data.fly.fv, 100*power_ratio.(animal_names{n}), ...
        'Color', cc(n,:), 'LineWidth', 2);
end
xlabel('Frequency (hz')
ylabel('Power factor (%)')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax, 'XScale', 'log', 'YScale', 'linear')
xlim([0.1 20])

ylim([-10^(2) 10^(6)])
symlog(ax,'y',1)

yline(0, 'k--');

leg = legend(h.power, animal_names);
set(leg, 'Box', 'off')

end
