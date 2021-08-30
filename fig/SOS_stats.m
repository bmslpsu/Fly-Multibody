function [] = SOS_stats()
%% SOS_stats:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL');

%% Make stats table
clearvars -except PATH FILE ALL
vI = 2;

clss = ["HeadFree", "BodyFixed"];
% trf = ["ref2head", "ref2head"];
trf = ["err2head", "err2head"];
prop = ["gain", "phase", "error", "IO_coherence"];
n_clss = length(clss);
n_prop = length(prop);

tbl_vars = ["condition", "fly", "frequency"];
clss_count = 0;
STATS = cell(n_clss,1);
for n = 1:n_clss
    IOFv = ALL.(clss(n)).FRF_data.IOFv{vI};
    IOFv = round(IOFv, 2);
    prop_map_all = [];
    for p = 1:n_prop
        prop_data = ALL.(clss(n)).FRF_data.(trf(n)).fly(vI).(prop(p));
        [n_freq, n_fly] = size(prop_data);
        prop_map = reshape(prop_data, [n_freq*n_fly 1]);
        if p == 1
            flyI = clss_count + (1:n_fly)';
            fly_map = cellfun(@(x) x*ones(n_freq,1), num2cell(flyI), 'UniformOutput', false);
            fly_map = cat(1, fly_map{:});
            IOFv_map = repmat(IOFv, [n_fly 1]);
            prop_map_all = [prop_map_all , [fly_map, IOFv_map , prop_map]];
        else
            prop_map_all = [prop_map_all ,  prop_map];
        end
    end
    T = splitvars(table(prop_map_all));
    cond = clss(n) + "_" + trf(n);
    cond_map = repmat(cond, [n_freq*n_fly 1]);
    T = [table(cond_map) , T];
    T.Properties.VariableNames = [tbl_vars ,  prop];
    STATS{n} = T;
    
    clss_count = clss_count + n_fly;
end
STATS_ALL = cat(1, STATS{:});

%% Frequency stats
freq = unique(STATS_ALL.frequency);
n_freq = length(freq);
P = [];
for f = 1:n_freq
    freqI = STATS_ALL.frequency == freq(f);
    T = STATS_ALL(freqI,:);
    for p = 1:n_prop
        P(f).(prop(p)) = anova1(T.(prop(p)), T.condition, 'off');
    end
end

end