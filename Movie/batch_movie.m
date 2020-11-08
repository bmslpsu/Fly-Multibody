function [] = batch_movie(rootdir, rootpat, vidFs)
%% batch_movie: runs a MakeMovie function for selected files
%
%   INPUT:
%       rootdir   	:   root directory
%       rootpat   	:   root pattern directory
%       vidFs       :   output video frame rate
%
%   OUTPUT:
%       -
%

export = true; % save the videos
% vidFs = 50;
% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250';
% rootpat = 'C:\Users\BC\Box\Git\Arena\Patterns';

[FILES, PATH] = uigetfile({'*.*'},'Select movie file', rootdir, 'MultiSelect','on');
FILES = cellstr(FILES);
[patfile, patpath] = uigetfile({'*.mat'},'Select pattern file', rootpat, 'MultiSelect','on');
pat = fullfile(patpath,patfile);

nfile = length(FILES);
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    main = fullfile(PATH, FILES{file});
    montage_SOS(main, pat, vidFs, export);
    %make_montage_passive_head(main, vidFs, export)
end
disp('ALL DONE')
end