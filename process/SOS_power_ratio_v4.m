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


P_body = minreal(MODEL.HeadFree.P_norm.body);
C_body = minreal(MODEL.HeadFree.C_norm.body);

P_head = minreal(MODEL.HeadFree.P_norm.head);
C_head = minreal(MODEL.HeadFree.C_norm.head);

% delay = 0.02;
% [sys.fly, data.fly, ~] = switch_system_v3(P_body, P_head, C_body, C_head, delay, delay, false);
[sys.fly, data.fly, ~] = plant_power(P_body, P_head, false);

%% Time constants
tau_body = 1 ./ P_body.Denominator{1}(2);
tau_head = 1 ./ P_head.Denominator{1}(2);
P_norm = tf(1, [tau_body 1]);

tau_all = linspace(0.01, tau_body, 200)';
[~,hI] = min(abs(tau_all - tau_head));
tau_all(hI) = tau_body;
tau_ratio = tau_body ./ tau_all;
% tau_ratio = round(tau_body./tau_all,5,'significant');
tau_leg = "tau = " + tau_ratio;
n_tau = length(tau_all);
for t = 1:n_tau
    name = ['tau_' num2str(t)];
    P_new = P_norm;
    P_new.denominator{1}(1) = tau_all(t);
    [sys_tau.(name), data_tau.(name)] = plant_power(P_norm, P_new, false); 
end

%% Time constant contour
clc
power_ratio_struct = structfun(@(x) x.power.ratio_1, data_tau, 'UniformOutput', false);
power_ratio = struct2cell(power_ratio_struct);
power_ratio = cat(2,power_ratio{:});

fv = data_tau.tau_1.fv;
[X,Y] = meshgrid(fv, tau_ratio);

fig = figure (201); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [4 3];
set(fig, 'Visible', 'on')
clear ax h
ax = subplot(1,1,1); hold on ; cla
% S = surf(X, 100*power_ratio', Y, 'EdgeColor', 'interp');
S = pcolor(X, 100*power_ratio', Y);
set(S, 'EdgeColor', 'interp')

colormap(cool)
cb = colorbar;
cb.Label.String = 'tau scale';
view(0,90)
rotate3d on

xlabel('Frequency (hz')
ylabel('Power increase (%)')
zlabel('tau scaled')
set(ax, 'Color', 'none', 'LineWidth', 0.75)
set(ax, 'XScale', 'log', 'YScale', 'linear')
xlim([0.1 20])
ylim([-50 1500])
% ax.YLim(1) = -50;

tau_dispI = round(logspace(0.1,2.3,7));
% tau_dispI(5) = hI;
tau_ratio_disp = round(tau_ratio(tau_dispI));
for n = tau_dispI
    h.tau(n) = plot(fv, 100*power_ratio_struct.(['tau_' num2str(n)]), 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
end

h.power(1) = plot(data.fly.fv, 100*data.fly.power.ratio_1, 'Color', [1 0 0], 'LineWidth', 2);
h.power(2) = plot(data.human.fv, 100*data.human.power.ratio_1, 'Color', 'g', 'LineWidth', 2);
uistack(h.power, 'top')

% leg_name = string(animal_names);
% leg = legend(h.power(1), 'fly');
leg = legend(h.power, 'fly', 'human');
set(leg, 'Box', 'off')

% freezeColors
% delete(S)
% delete(h.power)
% delete(h.tau)
% set(ax, 'YColor', 'w')

%% Animal #1: human head-eye
close all
clc

zeta = 0.5;
wn = 2*pi*0.3;
P_head = tf(wn^(2), [1 2*zeta*wn wn^(2)]);
% P_head = tf([wn^(2) 0], [1 2*zeta*wn wn^(2)]);
P_eye = tf(1, [0.13 1]);
[sys.human, data.human, ~] = plant_power(P_head, P_eye, true);

end
