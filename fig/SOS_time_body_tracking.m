function [] = SOS_time_body_tracking()
%% SOS_time_body_tracking:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select ellipse data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select registered data', root, 'MultiSelect','off');

datasets = ["Ellipse", "Registered"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'GRAND','FUNC','U','N');
end

%% Compare head free, head-fixed, body-fixed
clc
clearvars -except ALL FILE PATH

set_names = ["Ellipse", "Registered"];
state_names = ["body", "body"];
idx_state = [1 1 2];
cc = [0.9 0 0 ; 0 0.7 1];

n_set = length(state_names);
n_cond = ALL.Ellipse.N.vel;
vel = ALL.Ellipse.U.vel{1};

Fs = round(1 ./ mean(diff(ALL.Ellipse.GRAND.all(1).Time(:,:,1))));
time = ALL.Ellipse.GRAND.all(1).Time(:,:,1);
time_range = repmat( [0 20], [3 1]);

plot_fly = false;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 n_cond*1.2])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_cond,1);
for v = 1:n_cond
    data_span = [];
    [~,data_span(1)] = min(abs(time - time_range(v,1)));
    [~,data_span(2)] = min(abs(time - time_range(v,2)));
    data_span = data_span(1):data_span(2);
    
    time_span = (0:length(data_span)-1)' ./ Fs;
    func_time = ALL.Ellipse.FUNC{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ ALL.Ellipse.FUNC{v}.All.Fs;
    
  	ax(v) = subplot(n_cond,1,v); hold on
    title([num2str(vel(v)) ''], 'interpreter', 'none')
    
    axis tight
  	%set(ax(v,1), 'YTick', ALL.Ellipse.FUNC{v}.All.Amp*[0 1])
    %set(ax(v,1), 'YLim', ALL.Ellipse.FUNC{v}.All.Amp*[-1.3 1.3])
    plot(func_time_span, ALL.Ellipse.FUNC{v}.All.X(func_span), 'k', 'LineWidth', 0.5)
    
    for n = 1:n_set
        if plot_fly
           plot(time_span, ...
               squeeze(ALL.(set_names(n)).GRAND.fly_all(v).mean.State(data_span,idx_state(n),:)), ...
               'Color', [cc(n,:) 0.25], 'LineWidth', 0.25) 
        end
        
        [h.patch(v,n),h.line(v,n)] = PlotPatch(...
                ALL.(set_names(n)).GRAND.fly_stats(v).mean.State.mean(data_span,idx_state(n)),...
                ALL.(set_names(n)).GRAND.fly_stats(v).mean.State.std(data_span,idx_state(n)), ...
                time_span, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    uistack(h.line(v,:), 'top')
  	set(ax(v,:), 'XLim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
end
leg_label = state_names + "_" + set_names;
leg = legend(squeeze(h.line(1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
% leg.Position  = [0.16 0.38 0.63 0.1];
linkaxes(ax, 'xy')
set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
% set(ax, 'YLim', 15*[-1 1])
set(h.line, 'Marker', 'none','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.5)

% set(ax(1:end-1,:), 'XColor', 'none')
set(ax(1:end-1,:), 'XTick', [])

set(h.line(:,1), 'LineStyle', '--')
set(h.line, 'LineWidth', 1)

%% Calculate coherence between tracking types
fs = 100;
f = 0:0.05:50;
Cxy = cell(n_cond,1);
IOCxy = cell(n_cond,1);
IOFv = cell(n_cond,1);
for v = 1:n_cond
    X = squeeze(ALL.Ellipse.GRAND.all(v).State(:,1,:));
    Y = squeeze(ALL.Registered.GRAND.all(v).State(:,1,:));
    IOFv{v} = ALL.Ellipse.GRAND.all(v).IOFv(:,1,1);
    n_trial = size(X,2);
    for n = 1:n_trial
        cxy = mscohere(X(:,n), Y(:,n), [], [], f, fs);
        [~, IOcxy,~,~] = getfreqpeaks(f, cxy, [], IOFv{v}, 0.05, false);
        Cxy{v}(:,n) = cxy;
        IOCxy{v}(:,n) = IOcxy;
    end
end

%% Plot coherence
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 1.5 n_cond*1.2])
movegui(fig, 'center')
clear ax h

for v = 1:n_cond
    ax(v,1) = subplot(n_cond,1,v) ; cla ; hold on
        plot(IOFv{v}, IOCxy{v}, '.-', 'Color', [0.5 0.5 0.5 0.2], 'MarkerSize', 3)
        plot(IOFv{v}, mean(IOCxy{v},2), '.-', 'Color', 'k', 'MarkerSize', 8)
end
set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'YLim', [0 1], 'XScale', 'log', 'XLim', [0.2 14], 'XTick', [0.1 1 10])

end