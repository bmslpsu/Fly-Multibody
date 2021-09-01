function [] = SOS_fit_step()
%% SOS_fit_step:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL','MODEL');

%% Set initialization parameters
clearvars -except PATH FILE ALL MODEL
% sys = MODEL.HeadFree.G.head;
% sys

%% err2body & err2head with data head-free
clc
clear ax h H

H.body = MODEL.HeadFree.H.direct.body;
H.head = MODEL.HeadFree.H.direct.head;
H.gaze = H.head + H.body;
plot_names = fieldnames(H);
n_plot = length(plot_names);

cc = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1];

ts = 0.001;
fs = 1 / ts;
tt = (0:ts:1)';
u = zeros(size(tt));
I_start = round(0 / ts) + 1;
u(I_start:end) = 1;

% Remove high frequency components
fc = 25;
% [b, a] = butter(3, 30 / (fs/2), 'low');
lpf = tf(1, [(1 / (2*pi*fc)) 1]);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
fig.Position(3:4) = [4 3];
movegui(fig, 'center')
h = gobjects(n_plot,1);
ax(1) = subplot(1,1,1); cla ; hold on ; box on
xlabel('time (s)')
ylabel('°s^{-1}')
plot(tt, u, '--k', 'LineWidth', 1)
for n = 1:n_plot
    sys = H.(plot_names{n}) * lpf;
    y = lsim(sys, u, tt, 0);
    %y = filtfilt(b, a, y);
    h(n) = plot(tt, squeeze(y), 'Color', cc(n,:));
    
    if n == n_plot
        y_final = y(end);
        yline(y_final, '--', 'Color', 'g')
    end
end
uistack(h(end), 'top')
set(ax, 'Color', 'none', 'LineWidth', 0.75)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
set(ax, 'XLim', [-0.01 0.2])
% set(ax(2:end,:), 'XTick', [0.1 1 10])
set(ax, 'YLim', [-0.1 1.1])
set(h, 'LineWidth', 1)
set(ax(2:end-1,:), 'XTickLabel', [])
% set(ax(2:end,2:end), 'YTickLabel', [])

end