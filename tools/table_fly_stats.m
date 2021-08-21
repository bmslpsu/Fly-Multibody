function [stats] = table_fly_stats(DATA, group_clm, data_clms, time_avg)
%% table_fly_stats: 
%
%   INPUTS:
%       dataset    	: table containing singal_attributes objects
%       use_val     : group by this column of the table or set to false to not group
%       time_avg    : average in time (boolean)
%
%   OUTPUTS:
%       -
%
if nargin < 4
    time_avg = false;
    if nargin < 3
        group_clm = false;
    end
end

if length(group_clm) > 1 % vector externally given
    G = group_clm;
elseif ~group_clm % nio groups but from flies
    G = ones(size(DATA,1),1);
else % group by column in table
    G = DATA{:,group_clm};
end

[val_fly_group, val_group, fly_group] = findgroups(G, DATA.fly);
n_val = length(unique(val_group));
n_fly = length(unique(DATA.fly));

% Get data from table & store in structure
names = DATA.Properties.VariableNames(data_clms);
n_clm = length(data_clms);
for n = 1:n_clm
    if ~iscell(DATA.(names{n}))
        DATA.(names{n}) = num2cell(DATA.(names{n}));
    end
    stats.(names{n}) = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x, DATA.(names{n}), ...
        'UniformOutput', false), val_fly_group);
    stats.(names{n}) = splitapply(@(x) {cat(2,x{:})}, DATA.(names{n}), val_fly_group);
end

% Group by values in one table column
stats = structfun(@(x) splitapply(@(y) {y}, x, findgroups(val_group)), ...
    stats, 'UniformOutput', false);
stats = structfun(@(x) cat(2,x{:}), stats, 'UniformOutput', false);

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

% Mean in time (if specified)
if time_avg
    stats = structfun(@(x) cellfun(@(y) mean(y,1), x, 'UniformOutput', false), ...
        stats, 'UniformOutput', false);
end

% Fly stats
stats.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    stats, 'UniformOutput', false);

% Fly means
fnames = string(fieldnames(stats.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:n_val
        temp = {stats.fly_stats.(fnames(f))(:,v).mean};
        emptyI = cellfun(@(x) isempty(x), temp);
        temp = temp(~emptyI);
        stats.fly_mean.(fnames(f)){v} = cat(2, temp{:});
    end
end

% Grand stats
stats.val_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    stats.fly_mean, 'UniformOutput', false);

% Add combined data
stats.comb = comb;
end