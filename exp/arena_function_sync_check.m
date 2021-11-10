%% 
close all
clc

Fs = All.Fs;
T = All.time(end);
Fs_t = fix(Fs) / 1;
% Fs_t = Fs;
tintrp = (0:(1/Fs_t):T)';

[TRIG,PAT] = sync_pattern_trigger(t_p(1:end-1), data(1:end-1,2), T, ...
    data(1:end-1,1), false, [], false, false);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx+5);
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx+5);
pat_pos = interp1(pat_time, pat_pos, tintrp, 'linear');
pat_vel = central_diff(pat_pos, 1/Fs_t);

all_pos = interp1(All.time, All.X_step, tintrp, 'linear');
all_vel = central_diff(all_pos, 1/Fs_t);

% [Fv, all_mag_pos , all_phase_pos , ~] = FFT(tintrp, all_pos);
% [~, all_mag_vel , all_phase_vel , ~] = FFT(tintrp, all_vel);
% 
% [~, pat_mag_pos , pat_phase_pos , ~] = FFT(tintrp, pat_pos);
% [~, pat_mag_vel , pat_phase_vel , ~] = FFT(tintrp, pat_vel);

[~, Fv, all_mag_pos, all_phase_pos] = chirpz(all_pos, Fs_t, 0, Fs_t/2);
[~, ~, all_mag_vel, all_phase_vel] = chirpz(all_vel, Fs_t, 0, Fs_t/2);

[~, ~, pat_mag_pos, pat_phase_pos] = chirpz(pat_pos, Fs_t, 0, Fs_t/2);
[~, ~, pat_mag_vel, pat_phase_vel] = chirpz(pat_vel, Fs_t, 0, Fs_t/2);

fig = figure (100) ; clf
set(fig, 'Color', 'w')
ax(1) = subplot(4,1,1) ; hold on ; ylabel('Position (°)') ; xlabel('Time (s)')
    plot(All.time, All.X_step - All.X_step(1), 'k')
    plot(tintrp, pat_pos - pat_pos(1), 'r')
    xlim([0 T])

ax(2) = subplot(4,1,2) ; hold on ; ylabel('Position (°)')
    plot(Fv, all_mag_pos, 'k')
    plot(Fv, pat_mag_pos, 'r')
    plot(All.Freq, All.Amp, 'g.', 'MarkerSize', 12);

ax(3) = subplot(4,1,3) ; hold on ;  ylabel('Velocity (°/s)')
    plot(Fv, all_mag_vel, 'k')
    plot(Fv, pat_mag_vel, 'r')
    plot(All.Freq, 2*pi*All.Amp.*All.Freq, 'g.', 'MarkerSize', 12);
    
ax(4) = subplot(4,1,4) ; hold on ;  ylabel('Phase (°)')
    plot(Fv, rad2deg(all_phase_pos), 'k')
    plot(Fv, rad2deg(pat_phase_pos), 'r')
    plot(All.Freq, rad2deg(All.Phase), 'g.', 'MarkerSize', 12);
    
set(ax, 'LineWidth', 1, 'Box', 'off', 'YGrid', 'on')
set(ax(2:end), 'XLim', [0.05 max(All.Freq) + 2])
linkaxes(ax(2:end), 'x')

%% Replay
close all
clc

vI = 2;

func = replay.pos.body_error_panel(:,vI);
IOFv = replay.IOFreq{vI};
IOMag = replay.freq.pos.body_error_panel.IOmag(:,vI);

Fs = 400.63;
T = replay.time(end);
Fs_t = fix(Fs) / 1;
% Fs_t = Fs;
tintrp = (0:(1/Fs_t):T)';

[TRIG,PAT] = sync_pattern_trigger(t_p(1:end-1), data(1:end-1,2), T, ...
    data(1:end-1,1), false, [], false, false);

pat_time = PAT.time_sync(PAT.sync:PAT.end_idx+5);
pat_pos = 3.75*PAT.pos(PAT.sync:PAT.end_idx+5);
pat_pos = interp1(pat_time, pat_pos, tintrp, 'linear');
pat_vel = central_diff(pat_pos, 1/Fs_t);

all_pos = interp1(replay.time, func, tintrp, 'linear');
all_vel = central_diff(all_pos, 1/Fs_t);

[~, Fv, all_mag_pos, all_phase_pos] = chirpz(all_pos, Fs_t, 0, Fs_t/2);
[~, ~, all_mag_vel, all_phase_vel] = chirpz(all_vel, Fs_t, 0, Fs_t/2);

[~, ~, pat_mag_pos, pat_phase_pos] = chirpz(pat_pos, Fs_t, 0, Fs_t/2);
[~, ~, pat_mag_vel, pat_phase_vel] = chirpz(pat_vel, Fs_t, 0, Fs_t/2);

fig = figure (100) ; clf
set(fig, 'Color', 'w')
ax(1) = subplot(4,1,1) ; hold on ; ylabel('Position (°)') ; xlabel('Time (s)')
    plot(replay.time, func - func(1), 'k')
    plot(tintrp, pat_pos - pat_pos(1), 'r')
    xlim([0 T])

ax(2) = subplot(4,1,2) ; hold on ; ylabel('Position (°)')
    plot(Fv, all_mag_pos, 'k')
    plot(Fv, pat_mag_pos, 'r')
    plot(IOFv, IOMag, 'g.', 'MarkerSize', 12);

ax(3) = subplot(4,1,3) ; hold on ;  ylabel('Velocity (°/s)')
    plot(Fv, all_mag_vel, 'k')
    plot(Fv, pat_mag_vel, 'r')
    plot(IOFv, 2*pi*IOMag.*IOFv, 'g.', 'MarkerSize', 12);
    
ax(4) = subplot(4,1,4) ; hold on ;  ylabel('Phase (°)')
    plot(Fv, rad2deg(all_phase_pos), 'k')
    plot(Fv, rad2deg(pat_phase_pos), 'r')
    %plot(All.Freq, rad2deg(All.Phase), 'g.', 'MarkerSize', 12);
    
set(ax, 'LineWidth', 1, 'Box', 'off', 'YGrid', 'on')
set(ax(2:end), 'XLim', [0.05 max(IOFv) + 2])
linkaxes(ax(2:end), 'x')



