function [] = batch_headtrack_edge(root, playback, export)
%% batch_headtrack_edge: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%       showpoint 	:  	save tracking annotations video
%
%   OUTPUT:
%       -
%

% export = false;
playback = 0;
root = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_forRoll';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head_edge');
mkdir(headdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')

    obj = headtracker_area_v1(vidData);
    obj = play_tracking(obj, playback, t_v);
    
    yaw = obj.yaw;
    roll_idx = obj.roll_idx;
    pivot = obj.pivot;
 	save(fullfile(headdir,FILES{file}),'-v7.3', 'yaw', 'roll_idx', 'pivot', 't_v')                                                                     
end
disp('ALL DONE')
end