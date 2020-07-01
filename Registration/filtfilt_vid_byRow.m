function [filtvid] = filtfilt_vid_byRow(vid,n,Fc,Fs,debug)
%% filtfilt_vid: applies zero-phase filter to each pixel in a greyscale video in time
%
%   INPUT:
%       vid         : video data
%       n           : filter order
%       Fc          : cutoff-frequency [Hz]
%       Fs          : sampling-frequency [Hz]
%       debug      	: show comparison
%
%   OUTPUT:
%       filtvid  	: filtered video data
%

if nargin < 5
    debug = false;
end

A = squeeze(vid);
dim = size(A);
% A = reshape(A,[dim(1)*dim(2) dim(3)]);
A = double(reshape(A,[dim(1)*dim(2) dim(3)]));
[b,a] = butter(n, Fc / (Fs/2), 'low');
filtvid = zeros(size(A));
for r = 1:size(A,1)
    filtvid(r,:) = filtfilt(b, a, A(r,:));
end
filtvid = uint8(reshape(filtvid,dim));

if debug
    fig = figure(1) ; clf
    set(fig, 'Color', 'w')
    ax(1) = subplot(2,1,1); axis image
    ax(2) = subplot(2,1,2); axis image
    set(ax, 'Box', 'on')
    for f = 1:dim(3)
        subplot(2,1,1) ; cla
            imshow(vid(:,:,f))
        subplot(2,1,2) ; cla
            imshow(filtvid(:,:,f))
        pause(0.001)
    end
end   

end