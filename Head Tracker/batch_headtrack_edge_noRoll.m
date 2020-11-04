function [] = batch_headtrack_edge_noRoll(root, playback, export)
%% batch_headtrack_edge_noRoll: runs head tracker for user selected video files
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

playback = 0;
export = false;
root = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250\registered';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

moviedir = fullfile(PATH, 'movie_edge');
headdir = fullfile(PATH,'tracked_head_edge');
mkdir(headdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'regvid','t_v')

    obj = headtracker_area_v2(regvid);
    if export
        obj = play_tracking(obj, playback, t_v, moviedir, FILES{file});
    else
        obj = play_tracking(obj, playback, t_v);
    end
    
%     cla; hold on
%     plot(t_v, obj.yaw, 'k', 'LineWidth', 1)
%     pause
%     close all
    
    yaw = obj.yaw;
    %roll_idx = obj.roll_idx;
    pivot = obj.pivot;
 	save(fullfile(headdir,FILES{file}),'-v7.3', 'yaw', 'pivot', 't_v')                                                                     
end
disp('ALL DONE')
end