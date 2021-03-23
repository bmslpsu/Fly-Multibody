function [] = SS_time_free_vs_fixed()
%% SOS_frf_free_vs_fixed:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[HeadFree_file,HeadFree_path] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

% [HeadFixed_file,HeadFixed_path] = uigetfile({'*.mat'}, ...
%     'Select head fixed data', root, 'MultiSelect','off');

[BodyFixed_file,BodyFixed_path] = uigetfile({'*.mat'}, ...
    'Select body fixed data', root, 'MultiSelect','off');

HeadFree = load(fullfile(HeadFree_path,HeadFree_file),'FRF_data','FLY_mean','FUNC','U','N');
% HeadFixed = load(fullfile(HeadFixed_path,HeadFixed_file),'FRF_data','FLY_mean','FUNC','U','N');
BodyFixed = load(fullfile(BodyFixed_path,BodyFixed_file),'FRF_data','FLY_mean','FUNC','U','N');

%% Body, wings
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file
clc

T = ["body", "dwba", "dwba", "head"];
cc = [0.9 0 0; 0.2 0.8 1; 0.1 0.9 0.2; 0 0.2 0.8];

% Set time range of data based on # of cycles
time = HeadFree.FLY_mean.time;
IOFv = HeadFree.U.freq{1};
n_freq = length(IOFv);
Fs = round(1 ./ mean(diff(time)));
Ts = 1 / Fs;
start_cycle = (3:3+n_freq-1)';
start_time = (start_cycle ./ IOFv);
time_range = [start_time.*ones(n_freq,1) , start_time.*ones(n_freq,1) + (4./IOFv)];
time_range = Ts*round(time_range./Ts);
func = HeadFree.FUNC;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_freq,2);
for v = 1:n_freq
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = (data_span(1):data_span(2))';
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = func{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ func{v}.All.Fs;
    
    subplot(ceil(n_freq/1),1,v); hold on ; title([num2str(IOFv(v)) ' Hz'])
    yyaxis left
        ax(v,1) = gca; set(ax(v,1), 'YColor', cc(1,:))
        axis tight
        %set(ax(v,1), 'YTick', func{v}.All.Amp*[0 1])
        %set(ax(v,1), 'YLim', func{v}.All.Amp*[-1.3 1.3])
        set(ax(v,1), 'YTick', func{v}.All.norm_vel*[0 1])
        set(ax(v,1), 'YLim', func{v}.All.norm_vel*[-1.3 1.3])
        %plot(func_time_span, func{v}.All.X(func_span), '-k', 'LineWidth', 0.5)
        plot(func_time_span, func{v}.All.dX(func_span), '-k', 'LineWidth', 0.5)
        
        [h.patch(v,1),h.line(v,1)] = PlotPatch(HeadFree.FLY_mean.grand_mean(v).(T(1))(data_span),...
                  HeadFree.FLY_mean.grand_std(v).(T(1))(data_span), time_span, ...
                  0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 0.75);
        [h.patch(v,3),h.line(v,3)] = PlotPatch(BodyFixed.FLY_mean.grand_mean(v).(T(3))(data_span),...
                  BodyFixed.FLY_mean.grand_std(v).(T(3))(data_span), time_span, ...
                  0, 1, cc(3,:), 0.7*cc(3,:), 0.2, 0.75);
        %[h.patch(v,4),h.line(v,4)] = PlotPatch(HeadFree.FLY_mean.grand_mean(v).(T(4))(data_span),...
                  %HeadFree.FLY_mean.grand_std(v).(T(4))(data_span), time_span, ...
                  %0, 1, cc(3,:), 0.7*cc(3,:), 0.2, 0.75);
    
    yyaxis right
        ax(v,2) = gca; set(ax(v,2), 'YColor', cc(2,:))
        [h.patch(v,2),h.line(v,2)] = PlotPatch(-HeadFree.FLY_mean.grand_mean(v).(T(2))(data_span),...
                  HeadFree.FLY_mean.grand_std(v).(T(2))(data_span), time_span, ...
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


%% Gaze, head
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file
clc

T = ["body", "gaze", "head", "head"];
cc = [0.9 0 0; 0.5 0.3 1; 0 0.4 1 ; 0 0.8 0.2];

% Set time range of data based on # of cycles
time = HeadFree.FLY_mean.time;
IOFv = HeadFree.U.freq{1};
n_freq = length(IOFv);
Fs = round(1 ./ mean(diff(time)));
Ts = 1 / Fs;
start_cycle = (3:3+n_freq-1)';
start_time = (start_cycle ./ IOFv);
time_range = [start_time.*ones(n_freq,1) , start_time.*ones(n_freq,1) + (3./IOFv)];
time_range = Ts*round(time_range./Ts);
func = HeadFree.FUNC;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_freq,1);
for v = 1:n_freq
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = (data_span(1):data_span(2))';
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = func{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ func{v}.All.Fs;
    
    ax(v,1) = subplot(ceil(n_freq/1),1,v); hold on ; title([num2str(IOFv(v)) ' Hz'])
        axis tight
        set(ax(v,1), 'YTick', func{v}.All.Amp*[0 1])
        set(ax(v,1), 'YLim', func{v}.All.Amp*[-1.3 1.3])
        %set(ax(v,1), 'YTick', func{v}.All.norm_vel*[0 1])
        %set(ax(v,1), 'YLim', func{v}.All.norm_vel*[-1.3 1.3])
        plot(func_time_span, func{v}.All.X(func_span), '-k', 'LineWidth', 0.5)
        
        [h.patch(v,1),h.line(v,1)] = PlotPatch(HeadFree.FLY_mean.grand_mean(v).(T(1))(data_span),...
                  HeadFree.FLY_mean.grand_std(v).(T(1))(data_span), time_span, ...
                  0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 0.75);
        [h.patch(v,2),h.line(v,2)] = PlotPatch(HeadFree.FLY_mean.grand_mean(v).(T(2))(data_span),...
                  HeadFree.FLY_mean.grand_std(v).(T(1))(data_span), time_span, ...
                  0, 1, cc(2,:), 0.7*cc(2,:), 0.2, 0.75);
        [h.patch(v,3),h.line(v,3)] = PlotPatch(HeadFree.FLY_mean.grand_mean(v).(T(3))(data_span),...
                  HeadFree.FLY_mean.grand_std(v).(T(2))(data_span), time_span, ...
                  0, 1, cc(3,:), 0.7*cc(3,:), 0.2, 0.75);
        [h.patch(v,4),h.line(v,4)] = PlotPatch(BodyFixed.FLY_mean.grand_mean(v).(T(4))(data_span),...
                  BodyFixed.FLY_mean.grand_std(v).(T(3))(data_span), time_span, ...
                  0, 1, cc(4,:), 0.7*cc(4,:), 0.2, 0.75);

  	set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
 end
leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east', 'Orientation', 'horizontal');
leg.Position = [0.2 0.96 0.5 0.02];
% linkaxes(ax(:,1), 'y')
% set(ax(:,1), 'YLim', 1.3*[-1 1])
% set(ax(:,1), 'YTick', 1*[0 1])
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')

%% head all fixed
clearvars -except HeadFree HeadFixed BodyFixed HeadFree_file HeadFree_file HeadFixed_file BodyFixed_file
clc

T = ["head"];
cc = [0 0.8 0.2];

% Set time range of data based on # of cycles
time = HeadFree.FLY_mean.time;
IOFv = HeadFree.U.freq{1};
n_freq = length(IOFv);
Fs = round(1 ./ mean(diff(time)));
Ts = 1 / Fs;
start_cycle = (1:1+n_freq-1)';
start_time = (start_cycle ./ IOFv);
time_range = [start_time.*ones(n_freq,1) , start_time.*ones(n_freq,1) + (6./IOFv)];
time_range = Ts*round(time_range./Ts);
func = HeadFree.FUNC;

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 9])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_freq,1);
for v = 1:n_freq
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = (data_span(1):data_span(2))';
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = func{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ func{v}.All.Fs;
    
    ax(v,1) = subplot(ceil(n_freq/1),1,v); hold on ; title([num2str(IOFv(v)) ' Hz'])
        axis tight
        set(ax(v,1), 'YTick', func{v}.All.Amp*[0 1])
        set(ax(v,1), 'YLim', func{v}.All.Amp*[-1.3 1.3])
        %set(ax(v,1), 'YTick', func{v}.All.norm_vel*[0 1])
        %set(ax(v,1), 'YLim', func{v}.All.norm_vel*[-1.3 1.3])
        plot(func_time_span, func{v}.All.X(func_span), '-k', 'LineWidth', 0.5)
        plot(time_span, BodyFixed.FLY_mean.fly(v).(T(1))(data_span,:), 'Color', [0.7*cc(1,:) 0.5])
        [h.patch(v,1),h.line(v,1)] = PlotPatch(BodyFixed.FLY_mean.grand_mean(v).(T(1))(data_span),...
                  BodyFixed.FLY_mean.grand_std(v).(T(1))(data_span), time_span, ...
                  0, 1, cc(1,:), 0.7*cc(1,:), 0.2, 0.75);

  	set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
 end
leg = legend(h.line(end,:), T, 'Box', 'off', 'Location', 'east', 'Orientation', 'horizontal');
leg.Position = [0.2 0.96 0.5 0.02];
% linkaxes(ax(:,1), 'y')
% set(ax(:,1), 'YLim', 1.3*[-1 1])
% set(ax(:,1), 'YTick', 1*[0 1])
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XColor', 'none')


end