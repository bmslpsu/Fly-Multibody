function [] = power_savings_surf()
%% power_savings_surf:
% calculates the power saving of moving a plant with fatser dynamics (head, smaller time constant)
% throught the same trajetory of a slower plant (body, larger time constant)

root = 'Q:\OneDrive - PSU\The Pennsylvania State University\Mongeau, Jean-Michel - Head_body_coordination\power savings'; % root location
[FILE,PATH] = uigetfile({'*.mat'}, 'Select model file', root, 'MultiSelect','off'); % select model file
load(fullfile(PATH,FILE),'MODEL'); % load in fit models

%% Simulate power savings for head in place of body
clearvars -except PATH FILE ALL MODEL
clc
close all

P_body = minreal(MODEL.HeadFree.P.body);
P_head = minreal(MODEL.HeadFree.P.head);

[sys.fly, data.fly, ~] = plant_power(P_body, P_head, true); % simulate power savings for head in place of body

%% Simulate power savings for distinct time constant ratios
tau_body = 1 ./ P_body.Denominator{1}(2); % nominal time constant ratio (body)
P_norm = tf(1, [tau_body 1]); % nominal plant (body)

% tau_all = flipud(linspace(0.00001, tau_body, 1000)'); % make array of time constants up to nominal value
tau_all = flipud(logspace(-5, log10(tau_body), 30)'); % make array of time constants up to nominal value
tau_ratio = tau_body ./ tau_all; % time constant ratio

showplot = false; % show the power savings at each iteration
n_tau = length(tau_all);
for t = 1:n_tau % each time constant ratio
    name = ['tau_' num2str(t)];
    P_new = P_norm;
    P_new.denominator{1}(1) = tau_all(t);
    [sys_tau.(name), data_tau.(name)] = plant_power(P_norm, P_new, showplot);  % simulate power savings
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

%% Map time constant ratio to size ratio
% Jean, this is the most critical part, it determines if fly morophology
% can predict power savings
size_ratio = (tau_ratio).^(1/4);

%% Find the time constant ratio correspondng to the morphology of a few animals
animal = [];
animal.size_ratio.human = (19e-2) / (23e-3); % Lee or Heymsfield / Bekerman
animal.size_ratio.marmoset = (3.5e-2) / (11e-3); % Pandey / Solomon
animal.size_ratio.blowfly = (12.23e-3) / (3e-3) ; % Bunchu / Stoffolano
animal.size_ratio.wasp = ((10^(1.1)) * 1e-3) / (((10^(0.7))^(1/3)) * 1e-3) ; % O’Donnell / O’Donnell
% animal.size_ratio.human = 4.35 / (7.18e-3); % Heymsfield / Heymsfield

fnames = fieldnames(animal.size_ratio);
n_animal = length(fnames);
for n = 1:n_animal
    [~,fI] = min(abs(animal.size_ratio.(fnames{n}) - size_ratio));
    animal.power_ratio.(fnames{n}) = power_ratio_f0(fI);
    animal.power_ratio_all.(fnames{n}) = power_ratio(:,fI);
end

%% Fly experimental power savings
[~,fI_fly] = min(abs(data.fly.fv - f0));
power_ratio_fly = data.fly.power.ratio_1;
power_ratio_f0_fly = power_ratio_fly(fI_fly); % power savings for fly
fly_size_ratio_1 = MODEL.morph.body.M / MODEL.morph.head.M; % if we use the mass ratio
fly_size_ratio_2 = MODEL.morph.body.L / MODEL.morph.head.L; % if we use the length ratio

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

h.surf = surf(fv_grid, size_ratio_grid, 100*power_ratio');
% h.surf = surfl(fv_grid, size_ratio_grid, 100*power_ratio');
% h.surf = waterfall(fv_grid, size_ratio_grid, power_ratio');

h.surf.EdgeColor = 'none';
h.surf.FaceAlpha = 0.8;

xlabel('Frequency (hz)')
ylabel('Size ratio')
zlabel('Power savings (%)')

set(ax, 'XScale', 'log')

% xlim([0 20])
% zlim([0 300])
ax.YLim(1) = 1;

% view(-34, 37)
view(-23, 30)

colormap(parula)
colorbar

hold on
animal.color = hsv(n_animal);
for n = 1:n_animal
    animal_size_ratio_grid = repmat(animal.size_ratio.(fnames{n}), [length(fv), 1]);
    h.animal(n) = plot3(fv, animal_size_ratio_grid, 100*animal.power_ratio_all.(fnames{n}), ...
        'Color', animal.color(n,:));
end
% set(h.animal, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 15)
fly_size_ratio_grid = repmat(fly_size_ratio_2, [length(fv), 1]);
% fly_size_ratio_grid = repmat(2.5^(1/4), [length(fv), 1]);
h.fly.exp = plot3(fv, fly_size_ratio_grid, 100*power_ratio_fly, ...
    'k-', 'MarkerEdgeColor', 'none', 'MarkerSize', 20, 'LineWidth', 1);

leg = legend([h.animal , h.fly.exp], [fnames ; {'fly'}]);
leg.Location = 'north';
leg.Position = [0.2001    0.5552    0.1823    0.2240];

%% Power savings vs size ratio
fig = figure (100); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'power savings')
fig.Position(3:4) = [4 3];
movegui(fig, 'center')
set(fig, 'Visible', 'on')

clear ax h
ax = subplot(1,1,1); hold on ; cla ; title(['Frequency: ' num2str(f0) ' hz'])
h.curve = plot(size_ratio, 100*power_ratio_f0, '-', 'LineWidth', 1);
h.fly.exp = plot([fly_size_ratio_1 , fly_size_ratio_2], ...
    100*[power_ratio_f0_fly , power_ratio_f0_fly], ...
    'k.-', 'MarkerEdgeColor', 'none', 'MarkerSize', 20, 'LineWidth', 1);

animal.color = hsv(n_animal);
for n = 1:n_animal
    h.animal(n) = plot(animal.size_ratio.(fnames{n}), 100*animal.power_ratio.(fnames{n}), ...
        'Color', animal.color(n,:));
end
set(h.animal, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 15)
leg = legend([h.animal , h.fly.exp], [fnames ; {'fly'}]);
leg.Position = [0.5979    0.4119    0.2432    0.1267];

ax.YLim(1) = -0.05*max(100*power_ratio_f0);
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
