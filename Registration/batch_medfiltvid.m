function [] = batch_medfiltvid(root, K)
%% batch_medfiltvid: low-pass filter videos in time, save to new folder
%
%   INPUT:
%       vid         : video data
%       K           : filter kernel [rows,columns,frames]
%
%   OUTPUT:
%       VID         : filtered video data
%

[FILES, PATH] = uigetfile({'*.mat'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

filtdir = fullfile(PATH,'vid_medfilt');
mkdir(filtdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')
    
    filtvid = squeeze(vidData);    
    filtvid = 1.5*medfilt3(filtvid, K, 'replicate');
    
    save(fullfile(filtdir,FILES{file}), '-v7.3', 'filtvid', 't_v', 'K')
end
disp('ALL DONE')

end