function [] = batch_register(root, crop_xy, trf_set)
%% batch_register: register selected videos
%
%   INPUT:
%       root        :   root directory
%       crop_xy  	:   crop rectangle or boolean to manually crop
%       trf_set     :   prior affine2D cell array for eaach frame
%
%   OUTPUT:
%       -
%

if nargin < 3
    trf_set = [];
    if nargin < 2
       crop_xy = []; 
    end
end

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

regdir = fullfile(PATH,'registered');
mkdir(regdir)

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,char(FILES(file))),'vidData','t_v')

    if file == 1
        if length(crop_xy) == 1
            [~, crop_xy] = imcrop(vidData(:,:,:,1));
        end
    end

    [regvid,trf] = register_video_translation(vidData, crop_xy, trf_set);

    save(fullfile(regdir,FILES{file}),'-v7.3','regvid','trf','t_v','crop_xy')
end
disp('ALL DONE')
end