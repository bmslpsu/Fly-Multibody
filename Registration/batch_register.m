function [] = batch_register(root)
%% batch_register: register selected videos
%
%   INPUT:
%       root   	:   root directory
%
%   OUTPUT:
%       -
%

% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

regdir = fullfile(PATH,'registered');
mkdir(regdir)

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')

    [regvid,trf] = register_video(vidData);

    save(fullfile(regdir,FILES{file}),'-v7.3','regvid','trf','t_v')
    %save(fullfile(regdir,FILES{file}),'-v7.3','regvid','trf')
end
disp('ALL DONE')
end