function [] = Hist_kinematics()
%% Hist_kinematics: compare kinematics of multiple data sets
root = 'E:\DATA\Magno_Data\Multibody';
[FILES,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
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

%% Compare stats from multiple experiments of the same type
clearvars -except ALL n_file FILES labels
clc

stats = cell(n_file,1);
use_val = true;
time_avg = true;
for n = 1:n_file
    stats{n} = val_fly_stats(ALL{n}, use_val, time_avg);
end
n_val = size(stats{1}.lwing,2);

%% Combine experiments
names = ["swba", "lwing", "rwing"];
n_name = length(names);
for f = 1:n_name
    fly_stats.(names(f)) = cell(1,n_val);
    all_stats.(names(f)) = cell(1,n_val);
    for n = 1:n_file
        temp_fly = stats{n}.fly_stats.(names(f));
        temp_all = stats{n}.comb.(names(f));
        for v = 1:n_val
            val_temp_fly = [temp_fly(:,v).mean]';
            val_temp_all = cat(2,temp_all{:,v});
            
            fly_stats.(names(f)){n,v}(:,1) = val_temp_fly(:);
            all_stats.(names(f)){n,v}(:,1) = val_temp_all(:);
        end
    end
end

comb_stats.fly = [];
comb_stats.all = [];
for f = 1:n_name
    comb_stats.fly.(names(f)) = cell(1,n_val);
    comb_stats.all.(names(f)) = cell(1,n_val);
    for v = 1:n_val
        comb_stats.fly.(names(f)){v} = cat(1, fly_stats.(names(f)){:,v});
        comb_stats.all.(names(f)){v} = cat(1, all_stats.(names(f)){:,v});
    end
end

comb_stats.fly.G = cell(1,n_val);
comb_stats.all.G = cell(1,n_val);
for v = 1:n_val
    comb_stats.fly.G{v} = cellfun(@(x,y) y*ones(size(x)), fly_stats.(names(1))(:,v), ...
        num2cell(1:n_file)', 'UniformOutput', false); 
    comb_stats.fly.G{v} = cat(1, comb_stats.fly.G{v}{:});
    
    comb_stats.all.G{v} = cellfun(@(x,y) y*ones(size(x)), all_stats.(names(1))(:,v), ...
        num2cell(1:n_file)', 'UniformOutput', false); 
    comb_stats.all.G{v} = cat(1, comb_stats.all.G{v}{:});
end

%% Stats
P.fly = nan(n_name, n_val);
P.all = nan(n_name, n_val);
for f = 1:n_name
    for v = 1:n_val
        %P.fly(f,v) = anova1(comb_stats.fly.(names(f)){v}, comb_stats.fly.G{v}, 'off');
        P.fly(f,v) = kruskalwallis(comb_stats.fly.(names(f)){v}, comb_stats.fly.G{v}, 'off');
        P.all(f,v) = kruskalwallis(comb_stats.all.(names(f)){v}, comb_stats.all.G{v}, 'off');
    end
end
    
%% Fly mean comparison
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_val 2*n_name])
ax = gobjects(n_name,n_val);
pp = 1;
spread = 0.2;
cc = jet(n_file);
rng(1)
for f = 1:n_name
    jitter = rand(size(comb_stats.fly.G{1})) * spread - (spread/2);
    G = comb_stats.fly.G{1};
    for v = 1:n_val
        ax(f,v) = subplot(n_name,n_val,pp); cla ; hold on
        title(['p = ' num2str(P.fly(f,v))])
            %b = boxchart(comb_stats.fly.(names(f)){v}, 'GroupByColor', comb_stats.fly.G{v});
            %set(b, 'BoxFaceAlpha', 0.5, 'MarkerStyle', 'none', 'Notch', 'off')
            bx = boxplot(comb_stats.fly.(names(f)){v}, comb_stats.fly.G{v}, ...
            'Width', 0.5, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 0.5);
        
            h = get(bx(5,:),{'XData','YData'});
            for c = 1:size(h,1)
               patch(h{c,1},h{c,2}, cc(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            end

            set(findobj(ax(f,v),'tag','Median'), 'Color', 'k','LineWidth', 1);
            set(findobj(ax(f,v),'tag','Box'), 'Color', 'none');
            set(findobj(ax(f,v),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
            set(findobj(ax(f,v),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
            ax(f,v).Children = ax(f,v).Children([end 1:end-1]); 

            if v == 1
               ylabel([char(names(f)) '(°)']) 
            end
            
            plot(G + jitter, comb_stats.fly.(names(f)){v}, '.', 'Color', 'r')
        
        pp = pp + 1;
    end
end
leg = legend(labels, 'Box', 'off', 'Interpreter', 'none', 'Orientation', 'horizontal');
leg.Position = [0.33 0.04 0.45 0.05];

set(ax, 'Color', 'none', 'LineWidth', 1.5, 'XColor', 'none', 'Box', 'off')
for f = 1:n_name
    linkaxes(ax(f,:), 'xy')
end
linkaxes(ax(2:3,:), 'xy')
set(ax(:,2:end), 'YColor', 'none')

% set(ax(1,:), 'YLim', [50 130])
% set(ax(2,:), 'YLim', 0.5*[-1 1])

%% All data points comparison
fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_val 2*n_name])
ax = gobjects(n_name,n_val);
h = gobjects(n_name, n_val, n_file);
pp = 1;
h_bins = {30:0.5:150, 15:0.5:75, 15:0.5:75, -15:0.1:15};
cmap = jet(n_file+2);
for f = 1:n_name
    for v = 1:n_val
        ax(f,v) = subplot(n_name,n_val,pp); cla ; hold on
        title(['p = ' num2str(P.all(f,v))])
        if v == 1
           xlabel([char(names(f)) ' (°)']) 
        end
        for n = 1:n_file
          	h(f,v,n) = histogram(all_stats.(names(f)){n,v}, h_bins{f}, 'Normalization', 'probability', ...
                'EdgeColor', 'none', 'FaceColor', cmap(n + (n-1)*2,:));
%           	histogram(stats{n}.comb.val_all.(names(f)){v}, h_bins{f}, 'Normalization', 'probability', ...
%                 'EdgeColor', 'none')
        end
        
        for n = 1:n_file
            mm = median(all_stats.(names(f)){n,v});
            xline(mm, 'Color', cmap(n + (n-1)*2,:))
        end
        ax(f,v).YLim(1) = -0.001;
        
        pp = pp + 1;
    end
end
leg = legend(labels, 'Box', 'off', 'Interpreter', 'none');
leg.Position = [0.25 0.95 0.45 0.05];

set(ax, 'Color', 'none', 'LineWidth', 1.5, 'XColor', 'k')
for f = 1:n_name
    linkaxes(ax(f,:), 'xy')
end
linkaxes(ax(2:3,:), 'xy')
set(ax(:,2:end), 'YColor', 'none')

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', 'Probability')

%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
dataset_names = [];
for n = 1:n_file
    dataset_names = [dataset_names '_' char(labels(n))];
end
filename = ['Distributions' dataset_names '.mat'];
save(fullfile(savedir, filename), 'stats', 'all_stats', 'comb_stats', 'labels', 'P', '-v7.3')
end

%% Stats function
function [stats] = val_fly_stats(dataset, use_val, time_avg)
    if nargin < 3
        time_avg = false;
        if nargin < 2
            use_val = true;
        end
    end

    DATA = dataset.DATA;
    
    if use_val
        [val_fly_group, val_group, ~] = findgroups(DATA{:,3}, DATA.fly);
        n_val = dataset.N{1,3};
    else
        [val_fly_group, val_group] = findgroups(DATA.fly);
        n_val = 1;
    end
    
    % Group variables
    stats.lwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.lwing, ...
        'UniformOutput', false), val_fly_group);
    stats.rwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) -x.position, DATA.rwing, ...
        'UniformOutput', false), val_fly_group);
    stats.dwba = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.dwba, ...
        'UniformOutput', false), val_fly_group);
    stats.swba = cellfun(@(x,y) x + y, stats.lwing, stats.rwing, ...
            'UniformOutput', false);
    stats.dwba = cellfun(@(x) x - median(x,'all'), stats.dwba, 'UniformOutput', false);
    
    if ~isobject([DATA.head{:}]') % head-fixed
        stats.head = cellfun(@(x) 0*x, stats.dwba, 'UniformOutput', false);
    else
        stats.head = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.head, ...
            'UniformOutput', false), val_fly_group);
    end
    stats.head = cellfun(@(x) x - median(x,'all'), stats.head, 'UniformOutput', false);
    
    if ~isobject([DATA.body{:}]') % body-fixed
        stats.body = cellfun(@(x) 0*x, stats.dwba, 'UniformOutput', false);
    else
        stats.body = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position - median(x.position), ...
            DATA.body, 'UniformOutput', false), val_fly_group);
    end
    
    if use_val
        stats = structfun(@(x) splitapply(@(y) {y}, x, val_group), stats, 'UniformOutput', false);
        stats = structfun(@(x) cat(2,x{:}), stats, 'UniformOutput', false);
    end
    
    % Store raw data
    comb = stats;
    fnames = string(fieldnames(comb));
    n_field = length(fnames);
    for f = 1:n_field
        for v = 1:n_val
            comb.all.(fnames(f)) = cat(2, comb.(fnames(f)){:});
            comb.val_all.(fnames(f)){v} = cat(2, comb.(fnames(f)){:,v});
            for n = 1:dataset.N.fly
                comb.fly_all.(fnames(f)){n,1} = cat(2, comb.(fnames(f)){n,:});
            end
        end
    end
    
    % Mean in time (if specified)
    if time_avg
        stats = structfun(@(x) cellfun(@(y) mean(y,1), x, 'UniformOutput', false), ...
            stats, 'UniformOutput', false);
    end
    
    % Fly stats
    stats.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
        stats, 'UniformOutput', false);

    fnames = string(fieldnames(stats.fly_stats));
    n_field = length(fnames);
    for f = 1:n_field
        for v = 1:n_val
            stats.fly_mean.(fnames(f)){v} = cat(2, stats.fly_stats.(fnames(f))(:,v).mean);
        end
    end
    
    % Grand stats
    stats.val_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
        stats.fly_mean, 'UniformOutput', false);
    
    % Add combined data
    stats.comb = comb;
end
