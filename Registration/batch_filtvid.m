function [] = batch_filtvid(root, n, Fc)
%% batch_filtvid: low-pass filter videos in time, save to new folder
%
%   INPUT:
%       vid         : video data
%       n           : # frames to filter
%
%   OUTPUT:
%       VID         : filtered video data
%

% root = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_30_HeadFixed';

[FILES, PATH] = uigetfile({'*.mat'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

filtdir = fullfile(PATH,'filtvid');
mkdir(filtdir)
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')
    
    Fs = round( 1 /mean(diff(t_v)) );

    filtvid = filtfilt_vid(3.5*vidData, n, Fc, Fs);
   	for f = 1:size(filtvid,3)
       filtvid(:,:,f) = medfilt2((filtvid(:,:,f)),9*[1 1]);
    end
    % filtvid = isolate_wing_vid(filtvid, false);
    % filtvid = 255*uint8(filtvid);
    
    save(fullfile(filtdir,FILES{file}), '-v7.3', 'filtvid', 't_v', 'n', 'Fc')
end
disp('ALL DONE')

end