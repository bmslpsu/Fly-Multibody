function [] = batch_register(root, crop_xy)
%% batch_register: register selected videos
%
%   INPUT:
%       root   	:   root directory
%       crop    : 
%
%   OUTPUT:
%       -
%

if nargin < 2
   crop_xy = []; 
end

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

regdir = fullfile(PATH,'registered');
mkdir(regdir)

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')

    if file == 1
        if length(crop_xy) == 1
            [~, crop_xy] = imcrop(vidData(:,:,:,1));
        end
    end

    [regvid,trf] = register_video(vidData, crop_xy);

    save(fullfile(regdir,FILES{file}),'-v7.3','regvid','trf','t_v','crop_xy')
end
disp('ALL DONE')
end