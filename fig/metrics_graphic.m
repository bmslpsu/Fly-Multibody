
clear ; close all ; clc
t = 0:(1/100):1;
f = 2;
u = 1*sin(2*pi*f*t);
y = 0.6*sin(2*pi*f*t - deg2rad(45));
e = u - y;

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches')
ax = subplot(1,1,1) ; cla ; hold on
plot(t, u, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
plot(t, y, 'Color', 'k', 'LineWidth', 1)
plot(t, e, 'Color', 'r', 'LineWidth', 1)
