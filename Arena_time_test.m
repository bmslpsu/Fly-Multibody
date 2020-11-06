
close all
clc

Fs = 50;

[TRIG,PAT] = sync_pattern_trigger(t_p, data(:,2), 20, data(:,1), false, [], false, true);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx);
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx);
pat_pos = interp1(pat_time, pat_pos, (0:(1/Fs):20)', 'linear');
pat_time = (0:(1/Fs):20)';
pat_vel = diff(pat_pos) / (1/Fs);

[Fv, Mag , Phs , FREQ] = FFT(pat_time, pat_vel);

subplot(2,1,1) ; hold on
plot(All.time, All.X_step - All.X_step(1), 'k')
plot(pat_time, pat_pos - pat_pos(1), 'r')

subplot(2,1,2) ; hold on
plot(Fv, All.mag.dX, 'k')
plot(Fv, Mag, 'r')
xlim([0 10])