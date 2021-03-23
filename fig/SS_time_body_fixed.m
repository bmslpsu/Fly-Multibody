function [] = SS_time_body_fixed()
%% SS_time_body_fixed:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N')

%% Body, wings
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

pI = [1 2];
T = ["head", "wing"];
n_plot = length(pI);
cc = [0.9 0 0 ; 0 0.4 1 ; 0.5 0.3 1];

Fs = round(1 ./ mean(diff(GRAND.all(1).Time(:,:,1))));
Ts = 1 / Fs;
time = GRAND.all(1).Time(:,:,1);
start_cycle = (3:3+N.freq-1)';
start_time = (start_cycle ./ U.freq{1});
time_range = [start_time.*ones(N.freq,1) , start_time.*ones(N.freq,1) + (4./U.freq{1})];
time_range = Ts*round(time_range./Ts);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.freq,1);
for v = 1:N.freq
    %data_span = (find(time==time_range(v,1)) : find(time==time_range(v,2)))';
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = data_span(1):data_span(2);
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = FUNC{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ FUNC{v}.All.Fs;
    
    ax(v,1) = subplot(ceil(N.freq/1),1,v); hold on ; title([num2str(U.freq{1}(v)) ' Hz'])
    axis tight
  	set(ax(v,1), 'YTick', FUNC{v}.All.Amp*[0 1])
    set(ax(v,1), 'YLim', FUNC{v}.All.Amp*[-1.3 1.3])
    plot(func_time_span, FUNC{v}.All.X(func_span), 'k', 'LineWidth', 0.5)
    for n = 1:n_plot
        %plot(time_span, squeeze(GRAND.all(v).refState(span,1,:)), 'k', 'LineWidth', 1)
        %plot(time_span, squeeze(GRAND.all(v).State(span,pI(n),:)), 'LineWidth', 0.5)
        [h.patch(v,n),h.line(v,n)] = PlotPatch(GRAND.fly_stats(v).mean.State.mean(data_span,pI(n)),...
                  GRAND.fly_stats(v).std.State.mean(data_span,pI(n)), time_span, ...
                  0, 1, cc(n,:), 0.7*cc(n,:), 0.2, 0.75);
    end
    set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
end
leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east');
leg.Position = [0.05 0.95 0.25 0.04];

set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')

%% Head saturation 
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc

pI = [1];
T = ["head"];
n_plot = length(pI);
cc = [0 0.8 0.2];

Fs = round(1 ./ mean(diff(GRAND.all(1).Time(:,:,1))));
Ts = 1 / Fs;
time = GRAND.all(1).Time(:,:,1);
start_cycle = (3:3+N.freq-1)';
start_time = (start_cycle ./ U.freq{1});
time_range = [start_time.*ones(N.freq,1) , start_time.*ones(N.freq,1) + (4./U.freq{1})];
time_range = Ts*round(time_range./Ts);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.freq,1);
for v = 1:N.freq
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = data_span(1):data_span(2);
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = FUNC{1}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ FUNC{1}.All.Fs;
    
    ax(v,1) = subplot(ceil(N.freq/1),1,v); hold on ; title([num2str(U.freq{1}(v)) ' Hz'])
    axis tight
  	set(ax(v,1), 'YTick', FUNC{v}.All.Amp*[-1 0 1])
    set(ax(v,1), 'YLim', FUNC{v}.All.Amp*1.8*[-1 1])
    plot(func_time_span, FUNC{v}.All.X(func_span), 'k', 'LineWidth', 0.5)
    %yline(15,'--r')
    for n = 1:n_plot
        %plot(time_span, squeeze(GRAND.all(v).refState(span,1,:)), 'k', 'LineWidth', 1)
        %plot(time_span, squeeze(GRAND.all(v).State(data_span,pI(n),:)), 'Color', [cc(n,:) 0.5], 'LineWidth', 0.25)
        plot(time_span, squeeze(GRAND.fly_all(v).mean.State(data_span,pI(n),:)), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        [h.patch(v,n),h.line(v,n)] = PlotPatch(GRAND.fly_stats(v).mean.State.mean(data_span,pI(n)),...
                  GRAND.fly_stats(v).mean.State.std(data_span,pI(n)), time_span, ...
                  1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    ax(v,1).XLim(1) = -0.05*ax(v,1).XLim(2);
end
% leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east');
% leg.Position = [0.92 0.84 0.07 0.13];

set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')

end