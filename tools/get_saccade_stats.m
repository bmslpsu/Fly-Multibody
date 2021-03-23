function [stats, SACCADE] = get_saccade_stats(dataset, name, use_val, time_win, get_names)
%% get_saccade_stats: 
%
%   INPUTS:
%       dataset         : table containing singal_attributes objects
%       name            : name of state in saccade table
%       use_val         : group by this coulum of the table or set to false to not group
%       time_win       	: time window
%       get_names       : names of varibales to pull out during saccade window
%
%   OUTPUTS:
%       -
%
if nargin < 5
    get_names = [];
    if nargin < 4
        time_win = 0.2;
        if nargin < 3
            use_val = true;
            if nargin < 2
                name = 'body';
            end
        end
    end
end
scd_name = [name '_saccade'];
get_names = string(cellstr(get_names));

clss = ["position", "velocity", "acceleration"];
n_clss = length(clss);
name_clss = name + "_" + clss;

n_get = length(get_names);
get_name_clss = cell(n_get,1);
for m = 1:n_get
   get_name_clss{m} = get_names(m) + "_" + clss;
end

n_file = size(dataset.DATA,1);
scd_obj = dataset.DATA.(scd_name);
scd_time = (-time_win:scd_obj{1}.Ts:time_win)';

SACCADE = [dataset.D , table(nan(n_file,1))];
SACCADE.Properties.VariableNames(end) = {'rate'};
% start_clss = 1 + size(SACCADE,2);
% SACCADE = [SACCADE, splitvars(table(num2cell(zeros(n_file,n_clss))))];
% SACCADE.Properties.VariableNames(start_clss:start_clss+n_clss-1) = cellstr(clss);
for n = 1:n_file
    if (n == 1) || ~mod(n,10) || (n == n_file)
        disp(n)
    end
    
    SACCADE.rate(n) = scd_obj{n}.rate; % store the saccade rate
    if scd_obj{n}.count > 0 % if any saccades
        for f = 1:n_clss
            % Get saccades from main signal
            [scds,~,~,~] = getSaccade(scd_obj{n}, dataset.DATA.(name){n}.(clss(f)), ...
                time_win, time_win, true);
            if strcmp(clss(f),'position') && strcmp(name,'body')
                scds = cellfun(@(x) x - mean(x(1:round(length(x)/4))), scds, 'UniformOutput', false);
            end
            SACCADE.(name_clss(f)){n} = cat(2,scds{:});
            
            % Get saccades from other signals
            for m = 1:n_get
                [scds,~,~,~] = getSaccade(scd_obj{n}, dataset.DATA.(get_names(m)){n}.(clss{f}), ...
                    time_win, time_win, true);
                if strcmp(clss(f),'position') && strcmp(get_names{m},'body')
                    scds = cellfun(@(x) x - mean(x(1:round(length(x)/3))), scds, 'UniformOutput', false);
                end
                SACCADE.(get_name_clss{m}(f)){n} = cat(2,scds{:});
            end
        end
    else
        for f = 1:n_clss
            SACCADE.(name_clss(f)){n} = nan*scd_time;
            for m = 1:n_get
                SACCADE.(get_name_clss{m}(f)){n} = nan*scd_time;
            end
        end
    end
end
% empty_rmv = cellfun(@(x) isempty(x), SACCADE.(clss(1)));
% SACCADE = SACCADE(~empty_rmv,:);

% Get grouping variables
if use_val
    [val_fly_group, val_group, fly_group] = findgroups(dataset.DATA{:,use_val}, dataset.DATA.fly);
    n_val = length(unique(val_group));
else
    [val_fly_group, val_group] = findgroups(dataset.DATA.fly);
    n_val = 1;
end
n_fly = length(unique(dataset.DATA.fly));

% Group by fly & category
for f = 1:n_clss
    stats.(name_clss(f)) = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x, ...
        SACCADE.(name_clss(f)), 'UniformOutput', false), val_fly_group);
    for m = 1:n_get
        stats.(get_name_clss{m}(f)) = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x, ...
            SACCADE.(get_name_clss{m}(f)), 'UniformOutput', false), val_fly_group);
    end
end

% Group by category is set
if use_val
    stats = structfun(@(x) splitapply(@(y) {y}, x, findgroups(val_group)), stats, 'UniformOutput', false);
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
        for n = 1:n_fly
            comb.fly_all.(fnames(f)){n,1} = cat(2, comb.(fnames(f)){n,:});
        end
    end
end

% Fly stats
stats.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    stats, 'UniformOutput', false);

% Fly means
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
stats.time = scd_time;
end