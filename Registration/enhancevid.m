function [VID1] = enhancevid(vid,n)
%% enhancevid: applies median filter to greayscale video in time
%
%   INPUT:
%       vid         : video data
%       n           : # frames to filter
%
%   OUTPUT:
%       VID         : filtered video data
%
% imgaussfilt

clc
clear all
close all

%%
vid = squeeze(vidData);
[y,x,nframe] = size(vid);

% [VID1] = medfilt_video(vid,3);
[VID2] = filtfilt_vid(vid,2,30,200);

%% 
m = 1; % image exponent
powervid = double(vid) .^ m;
filtvid = powervid;
parfor f = 1:nframe
    fprintf('%i \n', f)
    filtvid(:,:,f) = imgaussfilt(powervid(:,:,f), 2); % gaussian blur
    % filtvid(:,:,f) = medfilt2(filtvid(:,:,f), [3 3]); % 2D median filter
    % filtvid(:,:,f) = imsharpen(filtvid(:,:,f)); % sharpen
end
filtvid = uint8(filtvid);

%%
fig = figure(1) ; clf
set(fig, 'Color', 'w')
ax(1) = subplot(1,3,1); axis image
ax(2) = subplot(1,3,2); axis image
ax(3) = subplot(1,3,3); axis image
for f = 1:nframe
    subplot(1,3,1) ; cla
        imshow(vid(:,:,f))
    subplot(1,3,2) ; cla
        imshow(VID2(:,:,f))
    subplot(1,3,3) ; cla
        imshow(filtvid(:,:,f))
	pause(0.001)
end

end