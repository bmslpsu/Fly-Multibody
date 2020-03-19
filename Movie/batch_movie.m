function [] = batch_movie(rootdir, rootpat, head_debug)
%% batch_movie: runs a MakeMovie function for selected files
%
%   INPUT:
%       rootdir   	:   root directory
%       rootpat   	:   root pattern directory
%
%   OUTPUT:
%       -
%

% rootdir = 'E:\Experiment_SOS_v1';
% rootpat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

bodydir = fullfile(PATH,'tracked_body');
mkdir bodydir
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')
    
   	close all
    
    [bAngles,imgstats,initframe] = bodytracker(vidData, playback, head_debug);

    save(fullfile(bodydir,FILES{file}),'-v7.3','bAngles','imgstats','initframe','t_v')
end
disp('ALL DONE')
end