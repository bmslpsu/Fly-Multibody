function [] = SOS_passive()
%% SOS_passive:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE{1},PATH{1}] = uigetfile({'*.mat'},'Select head-free data', root, 'MultiSelect','off');
[FILE{2},PATH{2}] = uigetfile({'*.mat'},'Select passive data', root, 'MultiSelect','off');

datasets = ["HeadFree", "Passive"];
for n = 1:length(FILE)
   ALL.(datasets(n)) = load(fullfile(PATH{n},FILE{n}),'GRAND','FUNC','Replay','U','N');
end

%% Time domain
clc
clearvars -except ALL FILE PATH

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
movegui(fig, 'center')
clear ax h

cc.body = [0.9 0 0];
cc.head = [ 0 0.4 1];
vI = 2;
ax(1) = subplot(2,1,1); cla ; hold on ; ylabel('Body (°)')
    % plot(ALL.HeadFree.GRAND.fly_stats(vI).mean.Time.mean(:,1), ...
    %     ALL.HeadFree.GRAND.fly_stats(vI).mean.State.mean(:,1), ...
    %     'r', 'LineWidth', 1)
    
    all_body = squeeze(ALL.Passive.GRAND.all.refState);
    all_body = all_body - mean(all_body(1,:)) + ALL.Passive.Replay{1}.pos.body(1,vI);
    plot(squeeze(ALL.Passive.GRAND.all.Time), all_body, ...
        'Color', [cc.body 0.2], 'LineWidth', 0.25)
    plot(ALL.Passive.Replay{1}.time, ALL.Passive.Replay{1}.vel.body(:,vI), '--k', 'LineWidth', 0.5)
    
ax(2) = subplot(2,1,2); cla ; hold on ; ylabel('Head (°)')
    plot(ALL.HeadFree.GRAND.fly_stats(2).mean.Time.mean(:,1), ...
        ALL.HeadFree.GRAND.fly_stats(2).mean.State.mean(:,2), ...
        'k', 'LineWidth', 0.5)
    plot(squeeze(ALL.Passive.GRAND.all.Time), squeeze(ALL.Passive.GRAND.all.State), ...
        'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
    plot(ALL.Passive.GRAND.fly_stats(1).mean.Time.mean(:,1), ...
        ALL.Passive.GRAND.fly_stats(1).mean.State.mean(:,1), ...
        'Color', cc.head, 'LineWidth', 0.5)
    
    xlabel('Time (s)')

linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XLim', [-0.5 20])
set(ax(1), 'XColor', 'none')


%% Frequency domain
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 (2.5/2)*3])
movegui(fig, 'center')
clear ax h

cc.body = [0.9 0 0];
cc.head = [ 0 0.4 1];
cc.cohr = mean([cc.body ; cc.head],1);
vI = 2;
ax(1) = subplot(3,1,1); cla ; hold on ; ylabel('Body (°/s)')
    % plot(ALL.HeadFree.GRAND.fly_stats(vI).mean.Time.mean(:,1), ...
    %     ALL.HeadFree.GRAND.fly_stats(vI).mean.State.mean(:,1), ...
    %     'r', 'LineWidth', 1)
    
    plot(squeeze(ALL.Passive.GRAND.all.Fv), squeeze(ALL.Passive.GRAND.all.refMag), ...
        'Color', [cc.body 0.2], 'LineWidth', 0.25)
    plot(ALL.Passive.Replay{1}.Fv, ALL.Passive.Replay{1}.freq.vel.body.mag(:,vI), '--k', 'LineWidth', 0.5)
    
ax(2) = subplot(3,1,2); cla ; hold on ; ylabel('Head (°/s)')
    plot(ALL.HeadFree.GRAND.fly_stats(2).mean.Fv.mean(:,1), ...
        ALL.HeadFree.GRAND.fly_stats(2).mean.Mag.mean(:,2), ...
        'k', 'LineWidth', 0.5)
    plot(squeeze(ALL.Passive.GRAND.all.Fv), squeeze(ALL.Passive.GRAND.all.Mag), ...
        'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
    plot(ALL.Passive.GRAND.fly_stats(1).mean.Fv.mean(:,1), ...
        ALL.Passive.GRAND.fly_stats(1).mean.Mag.mean(:,1), ...
        'Color', cc.head, 'LineWidth', 0.5)
    ax(2).YLim(1) = -2;
    
    
ax(3) = subplot(3,1,3); cla ; hold on ; ylabel('Coherence')
    plot(squeeze(ALL.Passive.GRAND.all.Fv), squeeze(ALL.Passive.GRAND.all.Cohr), ...
        'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)    
        [h.patch, h.line] = PlotPatch(...
                ALL.Passive.GRAND.fly_stats(1).mean.Cohr.mean(2:end,1),...
                ALL.Passive.GRAND.fly_stats(1).mean.Cohr.std(2:end,1), ...
                ALL.Passive.GRAND.fly_stats(1).mean.Fv.mean(2:end,1), ...
                1, 1, cc.cohr, 0.7*cc.cohr, 0.2, 1);
    
    ax(3).YLim = [-0.05 1];
    xlabel('Frequency (Hz)')

linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XLim', [0.1 12], 'XScale', 'log', 'XTick', [0.1 1 10])
set(ax(1), 'XColor', 'none')
    
end