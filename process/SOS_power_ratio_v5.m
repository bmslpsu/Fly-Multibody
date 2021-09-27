function [] = SOS_power_ratio_v5()
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

P_body = minreal(MODEL.HeadFree.P.body);
C_body = minreal(MODEL.HeadFree.P.body);

P_head = minreal(MODEL.HeadFree.P.head);
C_head = minreal(MODEL.HeadFree.C.head);

% delay = 0.02;
% [sys.fly, data.fly, ~] = switch_system_v3(P_body, P_head, C_body, C_head, delay, delay, false);
[sys.fly, data.fly, ~] = plant_power(P_body, P_head, true);

%% Time constants
tau_body = 1 ./ P_body.Denominator{1}(2);
tau_head = 1 ./ P_head.Denominator{1}(2);
P_norm = tf(1, [tau_body 1]);

tau_all = flipud(linspace(0.001, tau_body, 300)');
% [~,hI] = min(abs(tau_all - tau_head));
% tau_all(hI) = tau_body;
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
f_cut_body = P_body.Denominator{1}(2) / (2*pi);
f0 = 5*f_cut_body;
[~,fI] = min(abs(fv - f0));
f0 = fv(fI);
power_ratio_f0 = power_ratio(fI, :);

%% Animal metric
size_ratio = (tau_ratio).^(1/2);

animal = [];
% animal.mass_ratio.human = 4.35 / (7.18e-3); % Heymsfield / Heymsfield
animal.mass_ratio.human = (19e-2) / (23e-3); % Lee or Heymsfield / Bekerman
animal.mass_ratio.moneky = (3.5e-2) / (11e-3); % Pandey / Solomon
animal.mass_ratio.blowfly = (12.23e-3) / (3e-3) ; % Bunchu / Stoffolano
animal.mass_ratio.wasp = ((10^(1.1)) * 1e-3) / (((10^(0.7))^(1/3)) * 1e-3) ; % O’Donnell / O’Donnell

fnames = fieldnames(animal.mass_ratio);
n_animal = length(fnames);
for n = 1:n_animal
    [~,fI] = min(abs(animal.mass_ratio.(fnames{n}) - size_ratio));
    animal.power_ratio.(fnames{n}) = power_ratio_f0(fI);
end

%% Fly experimental predicition
[~,fI_fly] = min(abs(data.fly.fv - f0));
f0_fly = data.fly.fv(fI_fly);
power_ratio_f0_fly = data.fly.power.ratio_1(fI_fly);
fly_size_ratio_1 = MODEL.morph.body.M / MODEL.morph.head.M;
fly_size_ratio_2 = MODEL.morph.body.L / MODEL.morph.head.L;

%% Power savings vs size ratio
fig = figure (309); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power')
fig.Position(3:4) = [4 3];
movegui(fig, 'center')
set(fig, 'Visible', 'on')

clear ax h
ax = subplot(1,1,1); hold on ; cla ; title(['Frequency: ' num2str(f0) ' hz'])
h.curve = plot(size_ratio, 100*power_ratio_f0, '-', 'LineWidth', 1);
h.fly.exp = plot([fly_size_ratio_1 , fly_size_ratio_2], ...
    100*[power_ratio_f0_fly , power_ratio_f0_fly], ...
    'r.-', 'MarkerEdgeColor', 'none', 'MarkerSize', 20, 'LineWidth', 1);

animal.color = hsv(n_animal);
for n = 1:n_animal
    h.animal(n) = plot(animal.mass_ratio.(fnames{n}), 100*animal.power_ratio.(fnames{n}), ...
        'Color', animal.color(n,:));
end
set(h.animal, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 15)
leg = legend(h.animal, fnames);
leg.Position = [0.5979    0.4119    0.2432    0.1267];

ax.YLim(1) = -0.05*(max(100*power_ratio_f0));
xlim([0 10^2])

xlabel('Size ratio')
ylabel('Power savings (%)')
zlabel('tau scaled')
set(ax, 'Color', 'none', 'LineWidth', 0.75)
set(ax, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
set(ax, 'XScale', 'log', 'YScale', 'linear')

set(ax, 'xminorgrid', 'on')
set(ax, 'GridLineStyle', '-')

end
