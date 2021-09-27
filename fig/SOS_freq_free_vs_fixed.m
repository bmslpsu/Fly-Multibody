function [] = SOS_freq_free_vs_fixed()
%% SOS_freq_free_vs_fixed:
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

n_set = length(state_names);
n_cond = ALL.HeadFree.N.vel;
vel = ALL.HeadFree.U.vel{1};
% amp = ALL.HeadFree.U.amp{1};
% cc = distinguishable_colors(n_set);
plot_fly = false;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 n_cond*1.2])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_cond,1);
for v = 1:n_cond    
  	ax(v) = subplot(n_cond,1,v); hold on
    title([num2str(vel(v)) ''], 'interpreter', 'none')
    IOFv = flipud( ALL.HeadFree.FUNC{v}.All.Freq );
    Fv = ALL.(set_names(1)).GRAND.fly_stats(v).mean.Fv.mean(:,1);
    plot(ALL.HeadFree.FUNC{v}.All.Fv, ALL.HeadFree.FUNC{v}.All.mag.dX, 'k', 'LineWidth', 0.5)
    
    set(ax(v), 'YTick', [0 vel(v)])
    ax(v).YLim(1) = -10;
    for n = 1:n_set
        if plot_fly
           plot(time_span, ...
               squeeze(ALL.(set_names(n)).GRAND.fly_all(v).mean.State(data_span,idx_state(n),:)), ...
               'Color', [cc(n,:) 0.25], 'LineWidth', 0.25) 
        end
        
%         if n == 1
%             mag = ALL.(set_names(n)).GRAND.fly_stats(v).mean.Mag.mean(:,idx_state(n));
%             mag(mag < 20) = 0;
%         else
%             mag = ALL.(set_names(n)).GRAND.fly_stats(v).mean.Mag.mean(:,idx_state(n));
%         end
        
        [h.patch(v,n),h.line(v,n)] = PlotPatch(...
                ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOMag.mean(:,idx_state(n)),...
                0*ALL.(set_names(n)).GRAND.fly_stats(v).mean.IOMag.std(:,idx_state(n)), ...
                IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
    end
    uistack(h.line(v,:), 'top')
  	%set(ax(v,:), Lim', [-0.05*ax(v,1).XLim(2) ax(v,1).XLim(2)+0.05*ax(v,1).XLim(2)])
end
leg_label = state_names + "_" + set_names;
leg = legend(squeeze(h.line(1,:)), leg_label, ...
    'Box', 'off', 'interpreter', 'none', 'Orientation', 'vertical');
% leg.Position  = [0.16 0.38 0.63 0.1];
linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XLim', [0.2 20], 'XScale', 'log')
% set(ax, 'YLim', 15*[-1 1])
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1.5)

% set(ax(1:end-1,:), 'XColor', 'none')
set(ax(1:end-1,:), 'XColor', 'none')

end