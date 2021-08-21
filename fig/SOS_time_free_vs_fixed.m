function [] = SOS_time_free_vs_fixed()
%% SOS_time_free_vs_fixed:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select head-free data', root, 'MultiSelect','off');
% [FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select head-fixed data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'}, 'Select body-fixed data', root, 'MultiSelect','off');

datasets = ["HeadFree", "BodyFixed"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'GRAND','FUNC','U','N');
end

%% Compare head free, head-fixed, body-fixed
clc
clearvars -except ALL FILE PATH

set_names = ["HeadFree", "BodyFixed", "HeadFree"];
state_names = ["body", "head", "head"];
idx_state = [1 1 2];
cc = [0.9 0 0 ; 0 0.8 0.2 ; 0 0.4 1];

% set_names = ["BodyFixed"];
% state_names = ["head"];
% idx_state = [1];
% cc = [0 0.8 0.2];

n_set = length(state_names);
n_cond = ALL.HeadFree.N.vel;
vel = ALL.HeadFree.U.vel{1};
% amp = ALL.HeadFree.U.amp{1};
% cc = distinguishable_colors(n_set);

Fs = round(1 ./ mean(diff(ALL.HeadFree.GRAND.all(1).Time(:,:,1))));
time = ALL.HeadFree.GRAND.all(1).Time(:,:,1);
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
    func_time = ALL.HeadFree.FUNC{v}.All.time;
    func_span = [];
    [~,func_span(1)] = min(abs(func_time - time_range(v,1)));
    [~,func_span(2)] = min(abs(func_time - time_range(v,2)));
    func_span = func_span(1):func_span(2);
    func_time_span = (0:length(func_span)-1)' ./ ALL.HeadFree.FUNC{v}.All.Fs;
    
  	ax(v) = subplot(n_cond,1,v); hold on
    title([num2str(vel(v)) ''], 'interpreter', 'none')
    
    axis tight
  	%set(ax(v,1), 'YTick', ALL.HeadFree.FUNC{v}.All.Amp*[0 1])
    %set(ax(v,1), 'YLim', ALL.HeadFree.FUNC{v}.All.Amp*[-1.3 1.3])
    plot(func_time_span, ALL.HeadFree.FUNC{v}.All.X(func_span), 'k', 'LineWidth', 0.5)
    
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

end