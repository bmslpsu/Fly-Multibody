function [] = batch_filtvid_magno(root, n, Fc)
%% batch_filtvid: low-pass filter videos in time, save to new folder
%
%   INPUT:
%       vid         : video data
%       n           : # frames to filter
%
%   OUTPUT:
%       VID         : filtered video data
%

% root = 'H:\EXPERIMENTS\RIGID\Experiment_Static_Wave';

[FILES, PATH] = uigetfile({'*.mat'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

filtdir = fullfile(PATH,'vid_filt');
mkdir(filtdir)

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'regvid','t_v')
    
    %Fs = round( 1 /mean(diff(t_v)) );
    Fs = 160;
    
    filtvid = filtfilt_vid(1*regvid, n, Fc, Fs);
%     for n = 1:size(filtvid,3)
%        filtvid(:,:,n) = medfilt2(5*filtvid(:,:,n),8*[1 1]);
%     end
    % filtvid = isolate_wing_vid(filtvid, false);
    
    %save(fullfile(filtdir,FILES{file}), '-v7.3', 'filtvid', 't_v', 'n', 'Fc')
    save(fullfile(filtdir,FILES{file}), '-v7.3', 'filtvid', 'n', 'Fc')
end
disp('ALL DONE')

end