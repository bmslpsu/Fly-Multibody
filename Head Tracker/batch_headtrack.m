function [] = batch_headtrack(root, npoints, playback, showpoint)
%% batch_headtrack: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       playback    :   playback rate
%
%   OUTPUT:
%       -
%

% showpoint = true;
% npoints = 10;
% playback = 5;
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS\registered';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head');

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'regvid','t_v')

    [hAngles,cPoint,validity,ROI,initframe,finalframe] = headtracker(regvid, npoints, playback, showpoint);
    
    save(fullfile(headdir,FILES{file}),'-v7.3','hAngles','cPoint','validity',...
                                            'ROI','initframe','finalframe','t_v')
    
    close all
end
disp('ALL DONE')
end