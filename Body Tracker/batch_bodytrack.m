function [] = batch_bodytrack(root, playback)
%% batch_bodytrack: runs body tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       playback    :   playback rate
%
%   OUTPUT:
%       -
%

% playback = 10;
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

bodydir = fullfile(PATH,'tracked_body');

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')

    [bAngles,imgstats,initframe] = bodytracker(vidData, playback);

    % save(fullfile(bodydir,FILES{file}),'-v7.3','bAngles','imgstats','initframe','t_v')
    
    close all
end
disp('ALL DONE')
end