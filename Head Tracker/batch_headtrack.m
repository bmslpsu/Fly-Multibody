function [] = batch_headtrack(root, npoints, center, playback, showpoint)
%% batch_headtrack: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       npoint      :   # of points for tracker
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%       showpoint 	:  	logical >>> true = debug mode
%
%   OUTPUT:
%       -
%

showpoint = true;
npoints = 1;
center = [];
playback = 20;
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_Ramp\registered';
% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS\registered';
% root = 'H:\EXPERIMENTS\RIGID\Experiment_Static_Wave';

% center = "H:\EXPERIMENTS\RIGID\Experiment_Sinusoid\18.75\Vid\tracked_head\fly_5_trial_10_freq_1.mat";
% center = load(center);
% center = center.cPoint;

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head');
mkdir(headdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'regvid','t_v')

    [hAngles,cPoint,validity,ROI,initframe,finalframe] = headtracker(regvid, npoints, center, ...
                                                                            playback, showpoint);
        
    figure
    imagesc(validity)
    
    pause
    
    close all
    
 	save(fullfile(headdir,FILES{file}),'-v7.3','hAngles','cPoint','validity',...
                                            'ROI','initframe','finalframe','t_v')
%  	save(fullfile(headdir,FILES{file}),'-v7.3','hAngles','cPoint','validity',...
%                                             'ROI','initframe','finalframe')
                                                                       
end
disp('ALL DONE')
end