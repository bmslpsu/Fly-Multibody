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

%% Compare stats from multiple experiments of the same type
clearvars -except ALL n_file FILES labels
clc

% stats = cell(n_file,1);
% use_val = false;
% time_avg = true;
% for n = 1:n_file
%     stats{n} = val_fly_stats(ALL{n}, use_val, time_avg);
% end
% n_val = size(stats{1}.lwing,2);

%% Detect head saccades
close all
clc

% Rigid saccade detection parameters
scd_rigid.thresh = [60, 1, 2, 0];
scd_rigid.true_thresh = 250;
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
scd_rigid.boundThresh = [0.2 60];

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
scd_magno.boundThresh = [0.2 60];

showplot = false;
HEAD_SACCADE.rigid = get_saccades(ALL{1}, scd_rigid, showplot);
HEAD_SACCADE.magno = get_saccades(ALL{2}, scd_magno, showplot);


%% Save
savedir = 'E:\DATA\Magno_Data\Multibody\processed';
dataset_names = [];
for n = 1:n_file
    dataset_names = [dataset_names '_' char(labels(n))];
end
filename = ['Distributions' dataset_names '.mat'];
save(fullfile(savedir, filename), 'stats', 'all_stats', 'comb_stats', 'labels', 'P', '-v7.3')
end

%% Get saccades function
function [SACCADE] = get_saccades(dataset, scd, showplot)
    time = cellfun(@(x) x.time, dataset.DATA.head, 'UniformOutput', false);
    pos_data = cellfun(@(x) x.position, dataset.DATA.head, 'UniformOutput', false);
    n_file = dataset.N.file;
    
    SACCADE = [dataset.D , table(nan(n_file,1))];
    SACCADE.Properties.VariableNames(end) = {'rate'};
    SACCADE = [SACCADE , splitvars(table(num2cell(zeros(n_file,1))))];
    SACCADE.Properties.VariableNames(end) = {'head'};
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

        if showplot
            figure (1)
            pause
            close all
        end
    end
end
