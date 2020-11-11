
close all
clc

Fs = 200;
T = 20;

[TRIG,PAT] = sync_pattern_trigger(t_p, data(:,2), 20, data(:,1), false, [], false, false);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx);
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx);
pat_pos = interp1(pat_time, pat_pos, (0:(1/Fs):T)', 'linear');
pat_time = (0:(1/Fs):T)';
pat_vel = diff(pat_pos) / (1/Fs);


subplot(2,1,1) ; hold on
plot(pat_time, 3.75*(func - mean(func)), 'k')
plot(pat_time, pat_pos - mean(pat_pos), 'r')
xlim([0 T])

subplot(2,1,2) ; hold on
[Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);
plot(Fv, Mag, 'k')
[Fv, Mag , ~ , ~] = FFT(pat_time, 3.75*func);
plot(Fv, Mag, 'r')
xlim([0.1 10])

% subplot(2,1,2) ; hold on
% [Fv, Mag , ~ , ~] = FFT(All.time, All.X_step);
% plot(Fv, Mag, 'k')
% 
% [Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);
% plot(Fv, Mag, 'r')
% xlim([0 10])