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

delay = 0.02;
[sys.fly, data.fly, ~] = switch_system_v3(P_body, P_head, C_body, C_head, delay, delay, true);


%% Animal #1: fly (from morphology)
clearvars -except PATH FILE ALL MODEL equations sys data body_tau head_pc body_static head_static

MODEL.morph.body.M = 0.84185066; % [mg] body mass
MODEL.morph.head.M = 0.08957730; % [mg] head mass
MODEL.morph.body.R = 2.22; % [mm] body radius
MODEL.morph.head.R = 0.42; % [mm] head radius
% MODEL.morph.body.Jzz = 0.56650720; % [mg*mm^2] body inertia
% MODEL.morph.head.Jzz = 0.00619410; % [mg*mm^2] head inertia

% head_mass = MODEL.morph.head.M; % [kg] Heymsfield
% head_radius = MODEL.morph.head.R; % [m] Rinkevich
% head_inertia = (2/5)*head_mass*head_radius^(2); % [kg*m^2]
head_inertia = MODEL.morph.head.Jzz; % [kg*m^2]

% body_mass = MODEL.morph.body.M;
% body_radius = MODEL.morph.body.R;
% body_inertia = (2/5)*body_mass*body_radius^(2); % [kg*m^2]
body_inertia = MODEL.morph.body.Jzz; % [kg*m^2]

inertia_ratio = head_inertia / body_inertia;
damp_ratio = head_static * inertia_ratio^(1/4);

P2_new = tf(damp_ratio, [inertia_ratio*body_tau 1]);
C2_new = minreal( sys.fly.G2 / P2_new );

[sys.fly2, data.fly2, ~] = switch_system_v3(sys.fly.P1, P2_new, sys.fly.C1, C2_new, 0, 0, true);


%% Animal comparison
close all ; clc
power_ratio = structfun(@(x) x.power.ratio_1, data, 'UniformOutput', false);
animal_names = fieldnames(power_ratio);
na = length(animal_names);

fig = figure (200);
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [3 2];
set(fig, 'Visible', 'on')
clear ax h

cc = hsv(na);
ax = subplot(1,1,1); hold on ; cla
yline(1, 'k--');
for n = 1:na
    h.power(n) = plot(data.fly.fv, 100*power_ratio.(animal_names{n}), ...
        'Color', cc(n,:), 'LineWidth', 2);
end
xlabel('Frequency (hz')
ylabel('Power factor (%)')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax, 'XScale', 'log')
xlim([0.1 100])
ax.YLim(1) = -100;

leg = legend(h.power, animal_names);
set(leg, 'Box', 'off')

end