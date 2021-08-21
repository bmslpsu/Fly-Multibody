function [] = Head_saccade()
%% Head_saccade: compare saccade rate and dynamics of data sets
warning('off', 'signal:findpeaks:largeMinPeakHeight')
root = 'E:\DATA\Magno_Data\Multibody';
[FILES,PATH] = uigetfile({'*.mat'}, 'Select data files', root, 'MultiSelect','on');
FILES = cellstr(FILES);
n_file = length(FILES);
ALL = cell(n_file,1);
labels = string(zeros(n_file,1));
for n = 1:n_file
    ALL{n} = load(fullfile(PATH,FILES{n}),'DATA','D','I','U','N');
    filedata = textscan(char(FILES{n}), '%s', 'delimiter', '_');
    filedata = filedata{1};
    ALL{n}.name = [filedata{1} '_' filedata{2}];
    labels(n) = ALL{n}.name;
end

%% Detect head saccades
clearvars -except ALL n_file FILES labels
close all
clc

% Rigid saccade detection parameters
scd_rigid.thresh = [50, 1, 1.5, 70];
% scd_rigid.thresh = [60, 1, 1, 0];
scd_rigid.true_thresh = 300;
scd_rigid.Fc_detect = [15 nan];
scd_rigid.Fc_ss = [nan nan];
scd_rigid.amp_cut = 4;
scd_rigid.dur_cut = 0.1;
scd_rigid.direction = 0;
scd_rigid.direction = 0;
scd_rigid.pks = [];
scd_rigid.sacd_length = nan;
scd_rigid.min_pkdist = 0.2;
scd_rigid.min_pkwidth = 0.02;
scd_rigid.min_pkprom = 50;
scd_rigid.min_pkthresh = 0;
scd_rigid.boundThresh = [0.25 50];

% Magno saccade detection parameters
scd_magno.thresh = [40, 1, 2, 0];
scd_magno.true_thresh = 200;
scd_magno.Fc_detect = [15 nan];
scd_magno.Fc_ss = [nan nan];
scd_magno.amp_cut = 4;
scd_magno.dur_cut = 0.1;
scd_magno.direction = 0;
scd_magno.direction = 0;
scd_magno.pks = [];
scd_magno.sacd_length = nan;
scd_magno.min_pkdist = 0.2;
scd_magno.min_pkwidth = 0.02;
scd_magno.min_pkprom = 50;
scd_magno.min_pkthresh = 0;
scd_magno.boundThresh = [0.25 50];

showplot = false;
HEAD_SACCADE.rigid = get_saccades(ALL{1}, scd_rigid, showplot); % RIGID: trial 29
HEAD_SACCADE.magno = get_saccades(ALL{2}, scd_magno, showplot);

% STATS = structfun(@(x) table_fly_stats(x, 4, [5 7 8 10:14], false), HEAD_SACCADE, 'UniformOutput', false); % SS
STATS = structfun(@(x) table_fly_stats(x, 3, [4 6 7 9:13], false), HEAD_SACCADE, 'UniformOutput', false); % SOS

%% Saccades
clss = 'rigid';
cc = [0 0.8 0.4];
tt = HEAD_SACCADE.rigid.time{2}(:,1);
n_val = size(STATS.rigid.position,2);
val = ALL{1}.U{1,3}{1};

fig = figure (12); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 1.7*n_val 4])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,n_val);
for n = 1:n_val
    ax(1,n) = subplot(2, n_val, n); hold on ; title(num2str(val(n)))
        plot(tt, STATS.(clss).comb.val_all.position{n}, 'Color', [0.5 0.5 0.5 0.2], 'lineWidth', 0.2)
        [h.std(1,n), h.mean(1,n)] = PlotPatch(STATS.(clss).val_stats.position(n).mean, ...
            STATS.(clss).val_stats.position(n).std, tt, 1, 1, cc, 0.7*cc, 0.3, 1);

    ax(2,n) = subplot(2, n_val, n + n_val); hold on
        plot(tt, STATS.(clss).comb.val_all.velocity{n}, 'Color', [0.5 0.5 0.5 0.2], 'lineWidth', 0.2)
        [h.std(2,n), h.mean(2,n)] = PlotPatch(STATS.(clss).val_stats.velocity(n).mean, ...
            STATS.(clss).val_stats.velocity(n).std, tt, 1, 1, cc, 0.7*cc, 0.3, 1);
end
set(ax, 'Color', 'none', 'LineWidth', 1)
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
set(ax(1,:), 'YLim', 25*[-1 1])
set(ax(2,:), 'YLim', [-500 1200])
set(ax, 'XLim', 0.1*[-1 1])
set(ax, 'XTick', -0.1:0.05:0.1)

set(ax(:,2:end), 'YColor', 'none')
set(ax(1,:), 'XColor', 'none')
set(ax(2,2:end), 'XTickLabel', [])

%% Saccades: rigid vs magno
clss = ["rigid", "magno"];
n_clss = length(clss);
valI = 2;
cc = [0 0.8 0.2 ; 0 0.4 1];

fig = figure (13); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 3])
movegui(fig, 'center')
clear ax h
ax = gobjects(2,1);
ax(1,1) = subplot(2, 1, 1); hold on ; title(num2str(val(valI)))
    for n = 1:n_clss
        plot(tt, STATS.(clss(n)).fly_mean.position{valI}, 'Color', [0.7*cc(n,:) 0.5], 'lineWidth', 0.5)
        [h.std(1,1), h.mean(1,1)] = PlotPatch(STATS.(clss(n)).val_stats.position(valI).median, ...
            STATS.(clss(n)).val_stats.position(valI).std, tt, 1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
    end
    ylabel('Displacement (°)')

ax(2,1) = subplot(2, 1, 2); hold on
    for n = 1:n_clss
        plot(tt, STATS.(clss(n)).fly_mean.velocity{valI}, 'Color', [0.7*cc(n,:) 0.5], 'lineWidth', 0.5)
        [h.std(2,1), h.mean(2,1)] = PlotPatch(STATS.(clss(n)).val_stats.velocity(valI).median, ...
            STATS.(clss(n)).val_stats.velocity(valI).std, tt, 1, 1, cc(n,:), 0.7*cc(n,:), 0.3, 1);
    end
    ylabel('Velocity (°s^{-1})')
    xlabel('Time (s)')
    
uistack(h.std(1,1), 'top')
uistack(h.std(2,1), 'top')
uistack(h.mean(1,1), 'top')
uistack(h.mean(2,1), 'top')

set(ax, 'Color', 'none', 'LineWidth', 1)
linkaxes(ax, 'x')
set(ax(1,:), 'YLim', 20*[-1 1])
set(ax(2,:), 'YLim', [-200 800])
set(ax, 'XLim', 0.1*[-1 1])
set(ax, 'XTick', -0.1:0.05:0.1)
set(ax(1), 'YTick', -20:5:20)
set(ax(1,:), 'XColor', 'none')
align_Ylabels(fig)

%% Saccade dynamics
% close all ; clc
clss = 'rigid';
stat_names = ["amplitude", "peak_vel", "duration", "start_pos", "end_pos", "rate"];
n_stat = length(stat_names);

fig = figure (14); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_stat 2])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_stat);
for a = 1:n_stat
    dyn = STATS.(clss).comb.val_all.(stat_names(a));
    dyn_all = cat(2, dyn{:})';
    G = cellfun(@(x,y) y*ones(1,size(x,2)), dyn, num2cell(val'), 'UniformOutput', false);
    G = cat(2, G{:})';
    STATS.(clss).dynamics.(stat_names(a)) = [dyn_all, G];

    ax(1,a) = subplot(1,n_stat, a); hold on ; ylabel(num2str(stat_names(a)), 'Interpreter', 'none') 
        bx = boxplot(dyn_all, G, ...
            'Width', 0.5, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 0.5, 'Colors', 'k');

        hbx = get(bx(5,:),{'XData','YData'});
        for c = 1:size(hbx,1)
           patch(hbx{c,1},hbx{c,2}, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        end    
end
set(ax, 'Color', 'none', 'LineWidth', 1, 'Box', 'off')

%% Combine rigid & magno while collapsing frequencies
valI = 2;
for a = 1:n_stat
    I = STATS.rigid.dynamics.(stat_names(a))(:,2) == val(valI); 
    R = STATS.rigid.dynamics.(stat_names(a))(I,1);
    
    I = STATS.magno.dynamics.(stat_names(a))(:,2) == val(valI); 
    M = STATS.magno.dynamics.(stat_names(a))(I,1);
    
    STATS.dynamics.(stat_names(a)) = [R ; M];
    STATS.dynamics.P.(stat_names(a)) = ranksum(R, M);
    
    R_fly = cat(2,STATS.rigid.fly_mean.(stat_names(a)){valI})';
    M_fly = cat(2,STATS.magno.fly_mean.(stat_names(a)){valI})';
 	STATS.dynamics.fly.(stat_names(a)) = [R_fly ; M_fly];
    STATS.dynamics.fly.P.(stat_names(a)) = ranksum(R_fly, M_fly);
    %[~,STATS.dynamics.fly.P.(stat_names(a))] = ttest2(R_fly, M_fly);
    
    STATS.dynamics.G.(stat_names(a)) = ...
        categorical([ones(size(R)) ; 2*ones(size(M))], [1 2], {'rigid', 'magno'});
    STATS.dynamics.fly.G.(stat_names(a)) = ...
        categorical([ones(size(R_fly)) ; 2*ones(size(M_fly))], [1 2], {'rigid', 'magno'});
end


%% Dynamics figure
fig = figure (15); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 0.75*[2 2 2*n_stat 1.5])
movegui(fig, 'center')
clear ax h
cc = [0 0.8 0.2 ; 0 0.4 1];
ax = gobjects(1,n_stat);
for a = 1:n_stat
    ax(1,a) = subplot(1,n_stat, a); hold on ; ylabel(num2str(stat_names(a)), 'Interpreter', 'none')
        bx = boxplot(STATS.dynamics.(stat_names(a)), STATS.dynamics.G.(stat_names(a)), ...
            'Width', 0.5, 'Symbol', '', 'Whisker', 1.5, 'OutlierSize', 0.5, 'Colors', 'k');

        hbx = get(bx(5,:),{'XData','YData'});
        for c = 1:size(hbx,1)
           patch(hbx{c,1},hbx{c,2}, cc(c,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        axt = ax(1,a);
        set(findobj(axt,'tag','Median'), 'Color', 'w', 'LineWidth', 1);
        set(findobj(axt,'tag','Box'), 'Color', 'none');
        set(findobj(axt,'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(axt,'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        axt.Children = axt.Children([end 1:end-1]);
end
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'Box', 'off')
set(ax(1:end), 'XColor', 'none')
set(ax(1), 'YLim', [0 30])
set(ax(2), 'YLim', [0 1200])
set(ax(3), 'YLim', [0 0.1])
set(ax(4:5), 'YLim', 22*[-1 1], 'YTick', -20:10:20)
set(ax(6), 'YLim', [0 1])

%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
fname = 'Saccade';
for n = 1:length(FILES)
    filedata = textscan(FILES{n}, '%s', 'delimiter', '._');
    temp_name = [];
    for k = 1:5
        temp_name = [temp_name '_' char(filedata{1}(k))];
    end
    fname = [fname temp_name];
end
root = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [fname '.mat']), 'HEAD_SACCADE', 'STATS');

end

%% Get saccades function
function [SACCADE] = get_saccades(dataset, scd, showplot)
    time = cellfun(@(x) x.time, dataset.DATA.head, 'UniformOutput', false);
    pos_data = cellfun(@(x) x.position, dataset.DATA.head, 'UniformOutput', false);
    n_file = dataset.N.file;
    
    scd_win = 0.2;
    int_win = [];
    norm = true;
    
    SACCADE = [dataset.D , table(nan(n_file,1))];
    %SACCADE = splitvars(SACCADE);
    eI = size(SACCADE,2);
    SACCADE.Properties.VariableNames(eI) = {'rate'};
    SACCADE = [SACCADE , splitvars(table(num2cell(zeros(n_file,8))))];
    SACCADE.Properties.VariableNames(eI+1:end) = ...
        {'head', 'position', 'velocity', 'time', 'amplitude', 'peak_vel', 'duration', 'start_pos'};
    for n = 1:n_file
        disp(n)
        % Detect saccades
        head_scd = saccade_v1(pos_data{n}, time{n}, scd.thresh, scd.true_thresh, ...
                        scd.Fc_detect, scd.Fc_ss, scd.amp_cut, scd.dur_cut, ...
                        scd.direction, scd.pks, scd.sacd_length, scd.min_pkdist, ...
                        scd.min_pkwidth, scd.min_pkprom, scd.min_pkthresh, ...
                        scd.boundThresh, showplot);
        
      	% Store saccades
        SACCADE.head{n} = head_scd;
        SACCADE.rate(n) = head_scd.rate;
        
        if SACCADE.head{n}.count > 0
            [scds,~,~,~] = getSaccade(head_scd, head_scd.position, scd_win, int_win, norm);
            SACCADE.position{n} = cat(2,scds{:});

            [scds,~,scd_time,~] = getSaccade(head_scd, head_scd.velocity, scd_win, int_win, norm);
            SACCADE.velocity{n} = cat(2,scds{:});
            SACCADE.time{n} = cat(2,scd_time{:});
            
            SACCADE.amplitude{n} = (head_scd.SACD.Direction.*head_scd.SACD.Amplitude)';
            SACCADE.peak_vel{n} = (head_scd.SACD.Direction.*head_scd.SACD.PeakVel)';
            SACCADE.duration{n} = head_scd.SACD.Duration';
            SACCADE.start_pos{n} = (head_scd.SACD.Direction.*head_scd.SACD.StartPos)';
            SACCADE.end_pos{n} = (head_scd.SACD.Direction.*head_scd.SACD.EndPos)';
        else
            SACCADE.position{n} = [];
            SACCADE.velocity{n} = [];
            SACCADE.time{n} = [];
            SACCADE.amplitude{n} = [];
            SACCADE.peak_vel{n} = [];
            SACCADE.duration{n} = [];
            SACCADE.start_pos{n} = [];
            SACCADE.end_pos{n} = [];
        end

        if showplot
            figure (1)
            pause
            close all
        end
    end
end
