function [] = SOS_make_replay()
%% SOS_make_replay:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'DATA','FUNC','GRAND','FLY','U','N')

%% Create body error replay panel function & figures
clearvars -except FILE FUNC DATA ALL GRAND FLY D I U N root
clc
replay = [];

replay.time = GRAND.all(1).Time(:,:,1);
Ts = mean(diff(replay.time));
Fs = round(1 / Ts);
[b, a] = butter(3, 15 / (Fs/2), 'low');

replay.Fv = GRAND.all(1).Fv(:,:,1);
replay.pos.ref = nan(length(replay.time), 3);
replay.pos.body = nan(length(replay.time), 3);
replay.pos.body_error = nan(length(replay.time), 3);
for v = 1:N.vel
    func_intrp = interp1(FUNC{v}.All.time, FUNC{v}.All.X, replay.time, 'pchip');
    
    replay.IOFreq{v} = GRAND.fly_stats(v).mean.IOFv.mean(:,1);
    replay.pos.ref(:,v) = GRAND.fly_stats(v).mean.refState.mean(:,1);
    replay.pos.body(:,v) = GRAND.fly_stats(v).mean.State.mean(:,1);
    replay.pos.body_error(:,v) = replay.pos.ref(:,v) - replay.pos.body(:,v);
    replay.pos.body_error(:,v) = filtfilt(b, a, replay.pos.body_error(:,v));
    
    replay.pos.body_error_smooth(:,v) = func_intrp - replay.pos.body(:,v);
    
    offset = replay.pos.body_error(:,v) + -min(replay.pos.body_error(:,v)) + 1;
    replay.panel(:,v) = round( 96 * offset / 360) + 1;
	replay.pos.body_error_panel(:,v) = 3.75 * replay.panel(:,v);
    
    if v == 1
        fnames = fieldnames(replay.pos);
    end
    for f = 1:length(fnames)
        replay.vel.(fnames{f})(:,v) = central_diff(replay.pos.(fnames{f})(:,v), Ts);
        
        [~, ~, replay.freq.pos.(fnames{f}).mag(:,v)] = chirpz(replay.pos.(fnames{f})(:,v), Fs, 0, Fs/2);
        [~, replay.freq.pos.(fnames{f}).IOmag(:,v)] = getfreqpeaks(replay.Fv, replay.freq.pos.(fnames{f}).mag(:,v), ...
            replay.freq.pos.(fnames{f}).mag(:,v), replay.IOFreq{v}, [], false);
        
        [~, ~, replay.freq.vel.(fnames{f}).mag(:,v)] = chirpz(replay.vel.(fnames{f})(:,v), Fs, 0, Fs/2);
        [~, replay.freq.vel.(fnames{f}).IOmag(:,v)] = getfreqpeaks(replay.Fv, replay.freq.vel.(fnames{f}).mag(:,v), ...
            replay.freq.vel.(fnames{f}).mag(:,v), replay.IOFreq{v}, [], false);
    end
end

%% Replay frequency domain
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 4 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,1);
for v = 1:N.vel
    ax(v,1) = subplot(N.vel,1,v); hold on ; title([ num2str(U.vel{1}(v)) '°/s'])
        plot(replay.IOFreq{v}, replay.freq.vel.body_error_panel.IOmag(:,v), '-g');
        
        h.mag(v,1) = plot(FUNC{v}.All.Fv, FUNC{v}.All.mag.dX, 'Color', 'k');
        %h.mag(v,1) = plot(replay.Fv, replay.freq.vel.ref.mag(:,v), 'k');
        h.mag(v,2) = plot(replay.Fv, replay.freq.vel.body.mag(:,v), 'r');
        h.mag(v,3) = plot(replay.Fv, replay.freq.vel.body_error.mag(:,v), 'c');
        
        h.IOmag(v,1) = plot(replay.IOFreq{v}, replay.freq.vel.ref.IOmag(:,v), 'k');
        h.IOmag(v,2) = plot(replay.IOFreq{v}, replay.freq.vel.body.IOmag(:,v), 'r');
        h.IOmag(v,3) = plot(replay.IOFreq{v}, replay.freq.vel.body_error.IOmag(:,v), 'c');
        
        ax(v).YLim(1) = -5; 
end
xlabel('Frequency (Hz)')
ylabel('Velocity (°/s)')
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.2 20], 'XTick', [0.1 1 10])
set(ax, 'XScale', 'log')
set(ax(1:end-1), 'XColor', 'none')
linkaxes(ax, 'x')

set(h.mag, 'LineWidth', 0.75)
set(h.IOmag, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 11, 'LineStyle', 'none')

leg = legend(h.mag(1,:), 'Perturbation', 'Body', 'Body error replay', 'Box', 'off');
leg.Position = [0.1521    0.8702    0.3745    0.1104];

%% Replay replay time domain
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 6 5])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,1);
for v = 1:N.vel
    ax(v,1) = subplot(N.vel,1,v); hold on ; title([ num2str(U.vel{1}(v)) '°/s'])
        %plot(replay.IOFreq{v}, replay.freq.vel.body_error_panel.IOmag(:,v), '-g');
        
        h.mag(v,1) = plot(FUNC{v}.All.time, FUNC{v}.All.X, 'Color', 'k');
        %h.mag(v,1) = plot(replay.time, replay.freq.vel.ref.mag(:,v), 'k');
        h.mag(v,2) = plot(replay.time, replay.pos.body(:,v), 'r');
        h.mag(v,3) = plot(replay.time, replay.pos.body_error(:,v), 'c');
        h.mag(v,4) = plot(replay.time, replay.pos.body_error_smooth(:,v), 'm');
        h.mag(v,5) = plot(replay.time, replay.pos.body_error_panel(:,v) ...
            - mean(replay.pos.body_error_panel(:,v)), 'g');
end
xlabel('Time (s)')
ylabel('Position (°)')
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [-0.5 20], ...
    'XTick', 0:2:20)
set(ax(1:end-1), 'XColor', 'none')
linkaxes(ax, 'x')

set(ax(2), 'YLim', 80*[-1 1])

set(h.mag, 'LineWidth', 0.75)

leg = legend(h.mag(1,:), 'Perturbation', 'Body', 'Body error replay', 'panel', ...
    'Box', 'off', 'Orientation', 'horizontal');
leg.Position = [0.1157    0.9086    0.7292    0.1104];

%% Save replay data
filedata = textscan(FILE, '%s', 'delimiter', '_');
dataset_name = [];
for n = 1:6
    dataset_name = [dataset_name '_' char(filedata{1}(n))];
end
fname = ['replay' dataset_name];

savedir = fullfile(root,'replay');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'replay', 'FUNC', 'U', 'N');

%% Make position functions
for v = 1:N.vel
   func = replay.panel(:,v); % at 100 Hz
   func = interp1(replay.time, func, FUNC{v}.All.time, 'nearest'); % interpolate to 400.63 Hz for arena
   fname = ['position_function_replay_SOS_vel_' num2str(U.vel{1}(v)) '.mat'];
   save(fullfile(savedir, fname), 'func')
end

end