function [regvid,trf] = register_video(vid, crop_xy)
%% register_video: register each frame of a video to the same reference
%
%   INPUT:
%       vid   	:   raw video
%
%   OUTPUT:
%       vid   	:   registered video
%       trf   	:   2D affine transformation for each frame
%

if nargin < 2
    crop_xy = [];
end

vid = squeeze(vid); % In case it's 4D
dim = size(vid); % video size

% Crop video
if isempty(crop_xy)
    regvid = vid;
elseif length(crop_xy) == 2
    x_range = round(crop_xy(1):dim(2)-crop_xy(1));
    y_range = round(crop_xy(2):dim(1)-crop_xy(2));
    regvid = vid(y_range,x_range,:);
elseif length(crop_xy) == 4 % rectangle
    init_frame = imcrop(vid(:,:,1), crop_xy);
    regvid = uint8(zeros(size(init_frame)));
    for n = 1:dim(3)
       regvid(:,:,n) = imcrop(vid(:,:,n), crop_xy);
    end
else
    error('"crop_xy" must be 2x1 or 4x1')
end

% Get intial orientation
v1 = regvid(:,:,1); % initial frame
v2 = medfilt2(imopen(v1,strel('disk',15)),[7 7]);
v3 = imbinarize(v2);
L = regionprops(v3,'Orientation'); % initial angle
refangle = L.Orientation;
[~,flip] = findheading(v1, 0);
if flip
    refangle = refangle - 180;
end

% Set up optimizer
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 1e-3;
optimizer.MaximumIterations = 150;
optimizer.MaximumStepLength = 0.03;
optimizer.MinimumStepLength = 0.0002;

trf     = cell(dim(3),1); % store 2D affine transformations here
fixed   = double(squeeze(regvid(:,:,1)));
sz      = imref2d(size(fixed));

% Register each frame with respect to the first frame
tic
for n = 1:dim(3)
    fprintf([int2str(n) '\n'])
    frame = double(regvid(:,:,n));
    if n == 1
        trf{n} = imregtform(frame, fixed, 'rigid', optimizer, metric);
    else
        for jj = n-1:-1:1
            if(isRigid(trf{jj}))
                break
            end
        end
        trf{n} = imregtform(frame, fixed, 'rigid', optimizer, metric,...
                            'InitialTransformation', trf{jj});
    end
    reg = imwarp(frame, trf{n}, 'OutputView', sz);
    regvid(:,:,n) = imrotate(reg, 90-refangle, 'crop');
    fixed = (fixed*n + reg)/(n+1);
end
regvid = fliplr(regvid);
disp('DONE')
toc
end