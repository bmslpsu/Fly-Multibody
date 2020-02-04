function [vid] = flipvid(vid,method)
%% flipud_vid: switches the names of two variables in a .mat data file
%
%   INPUT:
%       vidin     	:   input video data
%       method      :   direction to flip
%                           'ud'    - flip vertically
%                           'lr'    - flip horizontally
%                           'udlr'	- flip vertically & horizontally
%
%   OUTPUT:
%       vidout     	:   output video data
%

allmeth = ["ud","lr","udlr"];
whichmethod = find(strcmp(method,allmeth)==1);

if isempty(whichmethod)
    error('"%s" is not a valid method',method)
end

vid = squeeze(vid);
dim = size(vid);

for frame = 1:dim(3)
    switch whichmethod
        case 1
            vid(:,:,frame) = flipud(vid(:,:,frame));
        case 2
            vid(:,:,frame) = fliplr(vid(:,:,frame));
        case 3
            vid(:,:,frame) = fliplr(flipud(vid(:,:,frame)));
    end
end

end