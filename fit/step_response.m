
close all
clc

vel = 34.6667*3.75;

pat_pos = 3.75*round( 96*(data(:,2) ./ 10) );
pat_vel = diff(pat_pos);
sI = find(pat_vel > 0, 1, 'first');
pat_start = pat_pos(sI);
time_start = t_p(sI);

pat_vel = vel*ones(size(pat_pos));
pat_vel(1:sI) = 0;

tt = t_v;
fs = 1 / mean(diff(tt));
body_pos =  bAngles - bAngles(1);
head_pos = head_data.angle - head_data.angle(1);

[b, a] = butter(3, 12 / (fs/2), 'low');
body_pos = filtfilt(b, a, body_pos);
head_pos = filtfilt(b, a, head_pos);
% pat_pos = filtfilt(b, a, pat_pos);

body_vel = central_diff(body_pos, 1/fs);
head_vel = central_diff(head_pos, 1/fs);

body_vel = filtfilt(b, a, body_vel);
head_vel = filtfilt(b, a, head_vel);

clear ax
ax(1) = subplot(2,1,1); hold on
    plot(t_p, pat_pos, 'g')
    plot(time_start, pat_start, 'mo')
    plot(t_v, body_pos, 'r')
    plot(t_v, head_pos, 'b')
    
ax(2) = subplot(2,1,2); hold on
    plot(t_p, pat_vel, 'g')
    plot(t_v, body_vel, 'r')
    plot(t_v, head_vel, 'b')

set(ax, 'XLim', [0 1], 'LineWidth', 1)
linkaxes(ax, 'x')