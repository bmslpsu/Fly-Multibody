function [] = power_savings_surf_fitC()
%% power_savings_surf:
% calculates the power saving of moving a plant with fatser dynamics (head, smaller time constant)
% throught the same trajetory of a slower plant (body, larger time constant)

root = 'E:\DATA\Magno_Data\Multibody\Processed'; % root location
[FILE,PATH] = uigetfile({'*.mat'}, 'Select model file', root, 'MultiSelect','off'); % select model file
load(fullfile(PATH,FILE),'MODEL'); % load in fit models

%% Simulate power savings for head in place of body
clearvars -except PATH FILE ALL MODEL
clc
close all

y = 4; % tau mapping

P_body = minreal(MODEL.HeadFree.P.body);
P_head = minreal(MODEL.HeadFree.P.head);

% fv_all = (0:0.5:20)';
fv_all = [0.1 , 1:1:20]';
fv_all = (0.1:0.05:20)';
[sys.fly, data.fly, ~] = plant_power(P_body, P_head, true, fv_all); % simulate power savings for head in place of body
sys.fly.tau_ratio = (1 / P_body.Denominator{1}(2)) / (1 / P_head.Denominator{1}(2));
sys.fly.size_ratio = sys.fly.tau_ratio^(1/y);

%% Simulate power savings for distinct time constant ratios
tau_body = 1 ./ P_body.Denominator{1}(2); % nominal time constant ratio (body)
P_norm = tf(1, [tau_body 1]); % nominal plant (body)

% tau_all = flipud(linspace(0.00001, tau_body, 1000)'); % make array of time constants up to nominal value
tau_start = tau_body/(20^y);
tau_all = flipud(logspace(log10(tau_start), log10(tau_body), 100)'); % make array of time constants up to nominal value
tau_ratio = tau_body ./ tau_all; % time constant ratio

showplot = false; % show the power savings at each iteration
n_tau = length(tau_all);
for t = 1:n_tau % each time constant ratio
    name = ['tau_' num2str(t)];
    P_new = P_norm;
    P_new.denominator{1}(1) = tau_all(t);
    [sys_tau.(name), data_tau.(name)] = plant_power(P_norm, P_new, showplot, fv_all);  % simulate power savings
    if showplot && (t==1)
        ylim([0 20])
    end
end

% get the power savings for each time constant ratio at all frequencies
power_ratio_struct = structfun(@(x) x.power.ratio_1, data_tau, 'UniformOutput', false);
power_ratio = struct2cell(power_ratio_struct);
power_ratio = cat(2,power_ratio{:});

fv = data_tau.tau_1.fv; % frequency vector
f_cut_body = 1 / (tau_body*2*pi); % the cutoff freqyency of the nminal plant
f0 = 1*f_cut_body; % pick a frequency to display the power savings 
% (this is somewhat arbitrary, as long as the frequency is large enough to
% show power savings)

[~,fI] = min(abs(fv - f0)); % get the closest frequency in the frequency vector
f0 = fv(fI);  % get the closest frequency in the frequency vector
power_ratio_f0 = power_ratio(fI, :); % pull out power savings for the display frequenncy

%% Fly experimental power savings
[~,fI_fly] = min(abs(data.fly.fv - f0));
power_ratio_fly = data.fly.power.ratio_1;
power_ratio_f0_fly = power_ratio_fly(fI_fly); % power savings for fly
% fly_size_ratio_1 = MODEL.morph.body.M / MODEL.morph.head.M; % if we use the mass ratio
% fly_size_ratio_2 = MODEL.morph.body.L / MODEL.morph.head.L; % if we use the length ratio
fly_size_ratio = mean(MODEL.morph.body.L + MODEL.morph.body.R) / mean(MODEL.morph.head.L + MODEL.morph.head.R);
fly_size_ratio = (MODEL.morph.body.M / MODEL.morph.head.M)^(1/3);
fly_size_ratio = MODEL.morph.body.L / MODEL.morph.head.L;

%% Map time constant ratio to size ratio
% Jean, this is the most critical part, it determines if fly morophology
% can predict power savings
dS = fly_size_ratio - sys.fly.size_ratio;
size_ratio = (tau_ratio).^(1/y);

%% Find the time constant ratio correspondng to the morphology of a few animals
animal = [];
% animal.size_ratio.robot = (2.8) / (248e-3); % Iyer / Iyer
animal.size_ratio.robot = (2e-2) / (2.3e-3); % Iyer / Iyer
animal.size_ratio.human = (19e-2) / (23e-3); % Lee or Heymsfield / Bekerman
animal.size_ratio.marmoset = (3.5e-2) / (11e-3); % Pandey / Solomon
animal.size_ratio.mouse = ((21.9e-3) / (5.281 * 1e-3)) + 0.2 ; % Islam / Howland (Mus musculus)
animal.size_ratio.fruit_fly_theo = fly_size_ratio;
animal.size_ratio.blowfly = (12.23e-3) / (3e-3) ; % Bunchu / Stoffolano
animal.size_ratio.wasp = ((10^(1.1)) * 1e-3) / (((10^(0.7))^(1/3)) * 1e-3) ; % O’Donnell / O’Donnell
% animal.size_ratio.crow = (100e-3 - 58e-3) / (15.157e-3); % Ludwig / Howland
animal.size_ratio.crow = ((18*2.54)*1e-2) / (100e-3 - 58e-3); % / Ludwig
animal.size_ratio.zebrafish = (30e-3) / (2e-3); % Collery / Collery
animal.size_ratio.critical = 3;
% animal.size_ratio.fruit_fly_theo = fly_size_ratio;
% animal.size_ratio.fruit_fly_exp = sys.fly.size_ratio;
% animal.size_ratio.fruit_fly = fly_size_ratio;
% animal.size_ratio.fruit_fly = 1.7;

% animal.size_ratio.human = 4.35 / (7.18e-3); % Heymsfield / Heymsfield

fnames = fieldnames(animal.size_ratio);
n_animal = length(fnames);
for n = 1:n_animal
    [~,fI] = min(abs(animal.size_ratio.(fnames{n}) - size_ratio));
    animal.power_ratio.(fnames{n}) = power_ratio_f0(fI);
    animal.power_ratio_all.(fnames{n}) = power_ratio(:,fI);
end

%% Set up surface plot
clc
[fv_grid, size_ratio_grid] = meshgrid(fv, size_ratio);

fig = figure (100); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power savings')
fig.Position(3:4) = [6 4];
movegui(fig, 'center')
set(fig, 'Visible', 'on')
clear ax h
ax = subplot(1,1,1); cla

norm = max(power_ratio, [], 'all');
% h.surf = surfl(fv_grid, size_ratio_grid, 100*power_ratio');
h.surf = surf(fv_grid, size_ratio_grid, power_ratio' ./ norm);
% h.surf = waterfall(fv_grid, size_ratio_grid, power_ratio');

% h.light = light;
% h.light.Position = [1 0 1];
% h.light.Style = 'local';
% material shiny

% h.surf.EdgeColor = [0.5 0.5 0.5];
h.surf.EdgeColor = 'none';
% h.surf.EdgeAlpha = 0.5;
h.surf.FaceAlpha = 1;
% h.surf.FaceColor = [0 0 0.5];
% h.surf.FaceColor = 'b';
h.surf.AmbientStrength = 0.5;

xlabel('Frequency (hz)')
ylabel('Size ratio')
zlabel('Power savings (%)')

set(ax, 'XScale', 'log')
set(ax, 'YScale', 'log')

set(ax, 'XTick', 10.^[-1 0 1])
xlim([0.1 20])

% xlim([0 20])
% zlim([0 300])
ax.YLim(1) = 1;

% view(-34, 37)
view(-23, 30)
view(-48, 35)

cmap = copper(1000);
cmap = cmap(150:end,:);
colormap(cmap)
colorbar

hold on
% animal.color = hsv(n_animal);
bg = [0.3 0.2 0];
animal.color = distinguishable_colors(n_animal, bg);
for n = 1:n_animal
    animal_size_ratio_grid = repmat(animal.size_ratio.(fnames{n}), [length(fv), 1]);
    h.animal(n) = plot3(fv, animal_size_ratio_grid, animal.power_ratio_all.(fnames{n}) ./ norm, ...
        'Color', animal.color(n,:));
end
set(h.animal, 'Marker', 'none', 'MarkerFaceColor', 'none', 'MarkerSize', 15, 'LineWidth', 2)
set(h.animal(end), 'LineStyle', '--', 'Color', [0.5 0.5 0.5])
% fly_size_ratio_grid = repmat(fly_size_ratio, [length(fv), 1]);
% fly_size_ratio_grid = repmat(2.5^(1/4), [length(fv), 1]);
% h.fly.exp = plot3(fv, fly_size_ratio_grid, 100*power_ratio_fly, ...
%     'r-', 'MarkerEdgeColor', 'none', 'MarkerSize', 20, 'LineWidth', 2);

leg = legend([h.animal], [fnames ; {'fly'}], 'Interpreter', 'none');
leg.Location = 'north';
leg.Position = [0.2001    0.5552    0.1823    0.2240];

% delete(leg)
% delete(h.animal)
% delete(h.surf)
% 
% set(ax, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none', 'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off')
% delete(leg)

%% Legend colors
figure (10)
clf ; hold on
for n = 1:n_animal
    h.animal_color(n) = plot([0 1], [n n], 'Color', animal.color(n,:), 'LineWidth', 2);
end
ylim([0 n_animal+1])
% set(gca, 'Color', 'none')
% set(gcf, 'Color', 'none')

end
