function [vid] = flipvid(vid)
%% register_video: register each frame of a video to the same reference
%
%   INPUT:
%       vid   	:   raw video
%
%   OUTPUT:
%       vid   	:   registered video
%       trf   	:   2D affine transformation for each frame
%

for frame = 1:size(regvid,3)
   regvid(:,:,frame) = flipud(regvid(:,:,frame));
end

end