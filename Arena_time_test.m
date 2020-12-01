
close all
clc

Fs = All.Fs;
T = All.time(end);
% T = 10;

[TRIG,PAT] = sync_pattern_trigger(t_p(1:end-1), data(1:end-1,2), T, data(1:end-1,1), false, [], false, false);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx);
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx);
pat_pos = interp1(pat_time, pat_pos, (0:(1/Fs):T)', 'linear');
pat_time = (0:(1/Fs):T)';
pat_vel = diff(pat_pos) / (1/Fs);

% [Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);

subplot(2,1,1) ; hold on
plot(All.time, All.X_step - All.X_step(1), 'k')
plot(pat_time, pat_pos - pat_pos(1), 'r')
% plot(All.time, All.X_step - mean(All.X_step), 'k')
% plot(pat_time, pat_pos - mean(pat_pos), 'r')
xlim([0 T])

% subplot(2,1,2) ; hold on
% plot(All.Fv, All.mag.X_step, 'k')
% plot(Fv, Mag, 'r')
% xlim([0.1 10])

subplot(2,1,2) ; hold on
[Fv, Mag , ~ , ~] = FFT(All.time, All.X_step);
plot(Fv, Mag, 'k')

[Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);
plot(Fv, Mag, 'r')
xlim([0 15])