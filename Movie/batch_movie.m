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

[FILES, PATH] = uigetfile({'*.*'},'Select movie file', rootdir, 'MultiSelect','on');
[patfile, patpath] = uigetfile({'*.mat'},'Select pattern file', rootpat, 'MultiSelect','on');
pat = fullfile(patpath,patfile);

nfile = length(FILES);
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    main = fullfile(PATH, FILES{file});
    MagnoMontage_SOS_v2(main, pat, vidFs, export);
end
disp('ALL DONE')
end