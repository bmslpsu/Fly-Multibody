function [y_clean, x_clean, I_comb] = select_point_windows(y, x, xy_mode)
%% select_point_windows: 
% Maually select windows to remove from data vector, then remove them
%
%   INPUTS:
%       y           : y-axis data
%       x           : x-axis data
%       xy_mode   	: sets the mode - 'x', 'y',or 'xy'
%
%   OUTPUTS:
%       y_removed   : y-axis data with points removed
%       x_removed   : x-axis data with points removed
%       I_comb      : logical vector of removed points
%

if nargin < 3
   xy_mode = 'x'; 
end

I = nan(size(y));
I(:) = 1:length(y);
if isempty(x)
    x = I;
end

fig = figure;
set(fig, 'Color', 'w')
ax = subplot(1,1,1); hold on
plot(x, y, 'k', 'LineWidth', 1)

% Let user slect points
[x_win,y_win] = ginput;
checkN = mod(length(x_win),2);
if checkN
   warning('odd number of points selected, using last point as end point')
   x_win = [x_win ;  x_win(end)];
   y_win = [y_win ;  y_win(end)];
end
n_point = size(x_win,1);
n_window = n_point / 2;

% Reshape so xy pairs define windows
x_win = reshape(x_win, [2 n_window])';
y_win = reshape(y_win, [2 n_window])';

% Get the windows based on the mode
Iy_all = cell(n_window,1);
Ix_all = cell(n_window,1);
I_all = cell(n_window,1);
for n = 1:n_window
    Iy_all{n} = find((x >= min(x_win(n,:))) & (x <= max(x_win(n,:))));
    Ix_all{n} = find((y >= min(y_win(n,:))) & (y <= max(y_win(n,:))));

    switch xy_mode
        case 'xy'
            I_all{n} = intersect(Iy_all{n}, Ix_all{n});
            xp = [x_win(n,1) x_win(n,1) x_win(n,2) x_win(n,2)];
            yp = [y_win(n,1) y_win(n,2) y_win(n,2) y_win(n,1)];
        case 'x'
            I_all{n} = Iy_all{n};
            xp = [x_win(n,1) x_win(n,1) x_win(n,2) x_win(n,2)];
            yp = [ax.YLim(1) ax.YLim(2) ax.YLim(2) ax.YLim(1)];
        case 'y'
            I_all{n} = Ix_all{n};
            xp = [ax.XLim(1) ax.XLim(1) ax.XLim(2) ax.XLim(2)];
            yp = [y_win(n,1) y_win(n,2) y_win(n,2) y_win(n,1)];
    end
    patch(xp, yp, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
end

% Make logical vector to remove points
I_comb = cat(1, I_all{:});
I_log = any(I == I_comb', 2);

% Remove points
y_removed = y(I_log);
x_removed = x(I_log);
y_clean = y(~I_log);
x_clean = x(~I_log);

% Plot removed and cleanned data
plot(x_removed, y_removed, 'r.', 'MarkerSize', 8)
plot(x_clean, y_clean, 'g.', 'MarkerSize', 8)

end