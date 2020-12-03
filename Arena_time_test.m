
close all
clc

Fs = All.Fs;
T = All.time(end);
% T = 10;

[TRIG,PAT] = sync_pattern_trigger(t_p(1:end-1), data(1:end-1,2), T, data(1:end-1,1), false, [], false, false);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx+5) - 0*1/5000;
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx+5);
tinrtp = (0:(1/Fs):T)';
pat_pos = interp1(pat_time, pat_pos, tinrtp, 'linear');
pat_time = tinrtp;
% pat_vel = diff(pat_pos) / (1/Fs);
% pat_vel = [pat_vel(1) ; pat_vel];
pat_vel = central_diff(pat_pos, 1/Fs);
% [Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);

subplot(3,1,1) ; hold on
plot(All.time, All.X_step - All.X_step(1), 'k')
plot(pat_time, pat_pos - pat_pos(1), 'r')
% plot(All.time, All.X_step - mean(All.X_step), 'k')
% plot(pat_time, pat_pos - mean(pat_pos), 'r')
xlim([0 T])

% subplot(2,1,2) ; hold on
% plot(All.Fv, All.mag.X_step, 'k')
% plot(Fv, Mag, 'r')
% xlim([0.1 10])

subplot(3,1,2) ; hold on
[Fv, Mag , ~ , ~] = FFT(All.time, All.X_step);
plot(Fv, Mag, 'k')

[Fv, Mag , ~ , ~] = FFT(pat_time, pat_pos);
plot(Fv, Mag, 'r')
plot(All.Freq, All.Amp, 'g.', 'MarkerSize', 12);

xlim([0 20])
ylim([0 5])

subplot(3,1,3) ; hold on
[Fv, Mag , ~ , ~] = FFT(All.time, All.dX_step);
plot(Fv, Mag, 'k')

[Fv, Mag , ~ , ~] = FFT(pat_time, pat_vel);
plot(Fv, Mag, 'r')
plot(All.Freq, 2*pi*All.Amp.*All.Freq, 'g.', 'MarkerSize', 12);

xlim([0 20])
% ylim([0 5])