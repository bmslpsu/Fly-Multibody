function [stats] = val_fly_stats(DATA, use_val, time_avg)
%% val_fly_stats: 
%
%   INPUTS:
%       dataset    	: table containing singal_attributes objects
%       use_val     : group by this coulum of the table or set to false to not group
%       time_avg    : average in time (boolean)
%
%   OUTPUTS:
%       -
%
if nargin < 3
    time_avg = false;
    if nargin < 2
        use_val = true;
    end
end

if use_val
    [val_fly_group, val_group, fly_group] = findgroups(DATA{:,use_val}, DATA.fly);
    n_val = length(unique(val_group));
    %n_fly = length(unique(fly_group));
else
    [val_fly_group, val_group] = findgroups(DATA.fly);
    n_val = 1;
end
n_fly = length(unique(DATA.fly));

clss = ["position", "velocity", "acceleration"];

% Group variables
if ~isobject([DATA.dwba{:}]') % no wing data
    stats.lwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) 0*x.position, DATA.head, ...
        'UniformOutput', false), val_fly_group);
    stats.rwing = stats.lwing;
    stats.dwba = stats.lwing;
    stats.swba = stats.lwing;
else
    stats.lwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.lwing, ...
        'UniformOutput', false), val_fly_group);
    stats.rwing = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) -x.position, DATA.rwing, ...
        'UniformOutput', false), val_fly_group);
    stats.dwba = splitapply(@(x) {cat(2,x{:})}, cellfun(@(x) x.position, DATA.dwba, ...
        'UniformOutput', false), val_fly_group);
    stats.swba = cellfun(@(x,y) x + y, stats.lwing, stats.rwing, ...
            'UniformOutput', false);
    stats.dwba = cellfun(@(x) x - median(x,'all'), stats.dwba, 'UniformOutput', false);
end    

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
        stats.fly_mean.(fnames(f)){v} = cat(2, stats.fly_stats.(fnames(f))(:,v).mean);
    end
end

% Grand stats
stats.val_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    stats.fly_mean, 'UniformOutput', false);

% Add combined data
stats.comb = comb;
end