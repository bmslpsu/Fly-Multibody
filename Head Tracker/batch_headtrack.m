function [] = batch_headtrack(root, playback)
%% batch_headtrack: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       playback    :   playback rate
%
%   OUTPUT:
%       -
%

% playback = 10;
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS\registered';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head');

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')

    [bAngles,imgstats,initframe] = headtracker(vidData, playback);

    save(fullfile(headdir,FILES{file}),'-v7.3','bAngles','imgstats','initframe','t_v')
    
    close all
end
disp('ALL DONE')
end