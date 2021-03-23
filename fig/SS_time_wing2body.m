function [] = SS_time_wing2body()
%% SS_time:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'FUNC','DATA','GRAND','FLY','D','I','U','N');

%% Body, wings
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc
close all
% [stats] = val_fly_stats(DATA, 4, false);

close all
clc
pI = [1 5 10 11];
T = ["body", "dwba", "lwing", "rwing"];
cc = [1 0 0; 0 0.6 1 ; 0 0.8 0.1; 0.74 0 1];

time = DATA.reference{1}.time;
Fs = DATA.reference{1}.Fs;
Ts = 1 / Fs;
start_cycle = (1:1+N.freq-1)';
start_time = (start_cycle ./ U.freq{1});
time_range = [start_time.*ones(N.freq,1) , start_time.*ones(N.freq,1) + (4./U.freq{1})];
time_range = Ts*round(time_range./Ts);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.freq,2);
for v = 1:N.freq
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
    
    subplot(ceil(N.freq/1),1,v); hold on ; title([num2str(U.freq{1}(v)) ' Hz'])
    yyaxis left % body and reference
        ax(v,1) = gca; set(ax(v,1), 'YColor', cc(1,:))
        axis tight
        %set(ax(v,1), 'YTick', FUNC{v}.All.Amp*[0 1])
        %set(ax(v,1), 'YLim', FUNC{v}.All.Amp*[-1.3 1.3])
        set(ax(v,1), 'YLim', FUNC{v}.All.norm_vel*[-1.3 1.3])
        plot(func_time_span, FUNC{v}.All.dX(func_span) , '-k', 'LineWidth', 0.5)
        pos_state = GRAND.fly_stats(v).mean.State.mean(data_span,pI(1));
        temp = singal_attributes(pos_state, time, 12);
        [h.patch(v,1),h.line(v,1)] = PlotPatch(temp.position,...
                  GRAND.fly_stats(v).std.State.mean(data_span,pI(1)), time_span, ...
                  0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 0.75);
                  
    yyaxis right % wings
        ax(v,2) = gca; set(ax(v,2), 'YColor', cc(2,:))
        [h.patch(v,2),h.line(v,2)] = PlotPatch(-GRAND.fly_stats(v).mean.State.mean(data_span,pI(2)),...
                  GRAND.fly_stats(v).std.State.mean(data_span,pI(2)), time_span, ...
                  0, 1, cc(2,:), 0.7*cc(2,:), 0.2, 0.75);
%         [h.patch(v,3),h.line(v,3)] = PlotPatch(GRAND.fly_stats(v).mean.State.mean(data_span,pI(3)),...
%                   GRAND.fly_stats(v).std.State.mean(data_span,pI(3)), time_span, ...
%                   0, 1, cc(3,:), 0.7*cc(3,:), 0.2, 0.75);
%         [h.patch(v,4),h.line(v,4)] = PlotPatch(-GRAND.fly_stats(v).mean.State.mean(data_span,pI(4)),...
%                   GRAND.fly_stats(v).std.State.mean(data_span,pI(4)), time_span, ...
%                   0, 1, cc(4,:), 0.7*cc(4,:), 0.2, 0.75);
              
  	set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
 end
leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east', 'Orientation', 'horizontal');
leg.Position = [0.2 0.96 0.5 0.02];
linkaxes(ax(:,2), 'y')
set(ax(:,2), 'YLim', 1.3*[-1 1])
set(ax(:,2), 'YTick', 1*[0 1])
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')


%% Body, wings test
clearvars -except FUNC DATA ALL GRAND FLY D I U N root
clc
close all

close all
clc
pI = [1 5 10 11];
T = ["body", "dwba"];
cc = [1 0 0; 0 0.6 1];

time = DATA.reference{1}.time;
Fs = DATA.reference{1}.Fs;
Ts = 1 / Fs;
start_cycle = (3:3+N.freq-1)';
start_time = (start_cycle ./ U.freq{1});
time_range = [start_time.*ones(N.freq,1) , start_time.*ones(N.freq,1) + (4./U.freq{1})];
time_range = Ts*round(time_range./Ts);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.freq,2);
for v = 1:N.freq
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
    
    subplot(ceil(N.freq/1),1,v); hold on ; title([num2str(U.freq{1}(v)) ' Hz'])
    yyaxis left % body and reference
        ax(v,1) = gca; set(ax(v,1), 'YColor', cc(1,:))
        axis tight
        %set(ax(v,1), 'YTick', FUNC{v}.All.Amp*[0 1])
        %set(ax(v,1), 'YLim', FUNC{v}.All.Amp*[-1.3 1.3])
        %set(ax(v,1), 'YLim', FUNC{v}.All.norm_vel*[-1.3 1.3])
        %plot(func_time_span, FUNC{v}.All.norm_vel * FUNC{v}.All.X(func_span) ./ max(abs(FUNC{v}.All.X(func_span))), '-k', 'LineWidth', 0.5)
        %plot(func_time_span, FUNC{v}.All.dX(func_span), '--k', 'LineWidth', 0.5)
        %plot(func_time_span, FUNC{v}.All.Amp * FUNC{v}.All.dX(func_span) ./ FUNC{v}.All.norm_vel, '-m', 'LineWidth', 0.5)
        %plot(func_time_span, -FUNC{v}.All.X(func_span), '--k', 'LineWidth', 0.5)
        plot(func_time_span, FUNC{v}.All.X(func_span) , '-k', 'LineWidth', 0.5)
        pos_state = GRAND.fly_stats(v).mean.State.mean(data_span,pI(1));
        temp = singal_attributes(pos_state, time, 12);
        [h.patch(v,1),h.line(v,1)] = PlotPatch(temp.acceleration,...
                  GRAND.fly_stats(v).std.State.mean(data_span,pI(1)), time_span, ...
                  0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 0.75);
                  
    yyaxis right % wings
        ax(v,2) = gca; set(ax(v,2), 'YColor', cc(2,:))
        [h.patch(v,2),h.line(v,2)] = PlotPatch(GRAND.fly_stats(v).mean.State.mean(data_span,pI(2)),...
                  GRAND.fly_stats(v).std.State.mean(data_span,pI(2)), time_span, ...
                  0, 1, cc(2,:), 0.7*cc(2,:), 0.2, 0.75);
              
  	set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
 end
leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east', 'Orientation', 'horizontal');
leg.Position = [0.2 0.96 0.5 0.02];
linkaxes(ax(:,2), 'y')
set(ax(:,2), 'YLim', 1.3*[-1 1])
set(ax(:,2), 'YTick', 1*[0 1])
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')
end