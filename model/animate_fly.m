function [] = animate_fly(dim, Fs, body_ang, head_ang)
%% animate_fly:
%
%   INPUTS:
%       body_ang :
%
%   OUTPUTS:
%       -
%

fig = figure (8);
set(fig, 'Color', 'k', 'Units', 'inches', 'Position', [2 2 9 9])
movegui(fig, 'center')
ax = subplot(6,1,1:4); cla ; hold on
axis image
axis off
ax_scale = 0.8*(dim.a1 + dim.a2);
axis(2.5*[-1 1 -1 1])
set(ax, 'Color', 'k')

dynamic_fly = virtual_fly([0 0], dim, 'b', 'r', 1);
% dynamic_fly = draw(dynamic_fly, 0, 0, 60, 60);
for n = 1:length(body_ang)
    cla
    dynamic_fly = draw(dynamic_fly, body_ang(n), head_ang(n), 60, 60);
    pause(1/Fs)
end

end