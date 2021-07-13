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

%% Animal #1: human head-eye
close all
clc

zeta = 0.5;
wn = 2*pi*0.3;
% P_head = tf(wn^(2), [1 2*zeta*wn wn^(2)]);
P_head = tf([wn^(2) 0], [1 2*zeta*wn wn^(2)]);

P_eye = tf(1, [0.13 1]);

[sys.human, data.human, ~] = plant_power(P_head, P_eye, true);

%% Time constants
tau_norm = 0.1554;
P_norm = tf(1, [tau_norm 1]);

tau_all = linspace(0.01, 0.1554, 200)';
tau_ratio = round(tau_norm./tau_all,3,'significant');
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
fig.Position(3:4) = [5 4];
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

h.power(1) = plot(data.fly.fv, 100*data.fly.power.ratio_1, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
% h.power(2) = plot(data.human.fv, 100*data.human.power.ratio_1, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

uistack(h.power, 'top')

xlabel('Frequency (hz')
ylabel('Power factor (%)')
zlabel('tau scaled')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax, 'XScale', 'log', 'YScale', 'linear')
xlim([0.1 20])
ylim([-50 1200])
% ax.YLim(1) = -50;

% for n = 1:n_tau
%     plot(fv, 100*power_ratio_struct.(['tau_' num2str(n)]), 'k', 'LineWidth', 2);
% end


% leg_name = string(animal_names);
leg = legend(h.power(1), 'fly');
set(leg, 'Box', 'off')


end
