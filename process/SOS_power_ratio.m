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

P_body = MODEL.HeadFree.P.body;
C_body = MODEL.HeadFree.C.body;

P_head = MODEL.HeadFree.P.head;
C_head = MODEL.HeadFree.C.head;

% P_body = MODEL.HeadFree.P_norm.body;
% C_body = MODEL.HeadFree.C_norm.body;
% 
% P_head = MODEL.HeadFree.P_norm.head;
% C_head = MODEL.HeadFree.C_norm.head;

head_pc = C_head.Numerator{1}(1) / P_head.Numerator{1}(2);
body_tau = 1 / P_body.Denominator{1}(2);

delay = 0.02;
[sys.fly, data.fly, ~] = switch_system_v3(P_body, P_head, C_body, C_head, delay, delay, true);

%% Get approximate cutoff frequncies of closed-loop responses and compare to plants
[body_cl_fcut] = compute_cutoff(data.fly.fv, data.fly.gain.H1, data.fly.phase.H1, 'low');
[head_cl_fcut] = compute_cutoff(data.fly.fv, data.fly.gain.H2, data.fly.phase.H2, 'high', true);

body_plant_fcut = -pole(sys.fly.P1)/(2*pi);
head_plant_fcut = -pole(sys.fly.P2)/(2*pi);

body_trf = body_cl_fcut / body_plant_fcut;
head_trf = head_cl_fcut / head_plant_fcut;

%% Animal #1: human head-eye
clearvars -except PATH FILE ALL MODEL equations sys data body_tau head_pc

% I_spehere = 2/5*m*r^(2)
head_mass = 5; % [kg] Heymsfield
head_radius = 0.194; % [m] Rinkevich
head_inertia = (2/5)*head_mass*head_radius^(2); % [kg*m^2]

eye_mass = 7.2 * 10^(-3); % [kg] Heymsfield
eye_radius = 0.024521; % [m] Howland
eye_inertia = (2/5)*eye_mass*eye_radius^(2); % [kg*m^2]

% eye_mass = 1;
ratio = eye_inertia / head_inertia;

P1 = sys.fly.P1;
C1 = sys.fly.C1;
G1 = P1 * C1;

P2 = tf(1, [ratio*body_tau 1]);
C2 = tf([head_pc 0],1);

[sys.human1, data.human1, ~] = switch_system_v3(P1, P2, C1, C2, 0, 0, true);

%% Animal #2: human body-head
clearvars -except PATH FILE ALL MODEL equations sys data body_tau head_pc

% I_spehere = 2/5*m*r^(2)
head_mass = 5; % [kg] Heymsfield
head_radius = 0.194; % [m] Rinkevich
head_inertia = (2/5)*head_mass*head_radius^(2); % [kg*m^2]

body_mass = 90;
body_radius = 0.194;
body_inertia = (2/5)*body_mass*body_radius^(2); % [kg*m^2]

% eye_mass = 1;
ratio = head_inertia / body_inertia;

P1 = sys.fly.P1;
C1 = sys.fly.C1;

P2 = tf(1, [ratio*body_tau 1]);
C2 = tf([head_pc 0],1);

[sys.human2, data.human2, ~] = switch_system_v3(P1, P2, C1, C2, 0, 0, false);

%% Animal #3: bee body-head
clearvars -except PATH FILE ALL MODEL equations sys data body_tau head_pc

% I_spehere = 2/5*m*r^(2)
head_mass = 10e-3; % [kg] O’Donnell
head_radius = 10^(1/3) * 10^(-3); % [m] O'Donnell
head_inertia = (2/5)*head_mass*head_radius^(2); % [kg*m^2]

body_mass = (10^1.5)*10^(-3);
body_radius = 2e-3;
body_inertia = (2/5)*body_mass*body_radius^(2); % [kg*m^2]

% eye_mass = 1;
ratio = head_inertia / body_inertia;

P1 = sys.fly.P1;
C1 = sys.fly.C1;

P2 = tf(1, [ratio*body_tau 1]);
C2 = tf([head_pc 0],1);

[sys.bee, data.bee, ~] = switch_system_v3(P1, P2, C1, C2, 0, 0, false);

%% Animal #4: bird
clearvars -except PATH FILE ALL MODEL equations sys data body_tau head_pc

% I_spehere = 2/5*m*r^(2)
head_mass = 0.1; % [kg] Heymsfield
head_radius = 0.024521; % [m] Rinkevich
head_inertia = (2/5)*head_mass*head_radius^(2); % [kg*m^2]

eye_mass = 7.2 * 10^(-3); % [kg] Heymsfield
eye_radius = 0.024521; % [m] Howland
eye_inertia = (2/5)*eye_mass*eye_radius^(2); % [kg*m^2]

% eye_mass = 1;
ratio = eye_inertia / head_inertia;

P1 = sys.fly.P1;
C1 = sys.fly.C1;

P2 = tf(1, [ratio*body_tau 1]);
C2 = tf([head_pc 0],1);

[sys.bird, data.bird, ~] = switch_system_v3(P1, P2, C1, C2, 0, 0, false);

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
yline(1, 'k--')
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