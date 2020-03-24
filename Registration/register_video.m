function [vid,trf] = register_video(vid)
%% register_video: register each frame of a video to the same reference
%
%   INPUT:
%       vid   	:   raw video
%
%   OUTPUT:
%       vid   	:   registered video
%       trf   	:   2D affine transformation for each frame
%

% Set up optimizer
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 1e-3;
optimizer.MaximumIterations = 150;
optimizer.MaximumStepLength = 0.03;
optimizer.MinimumStepLength = 0.0002;

vid = flipup(squeeze(vid)); % In case it's 4D
dim = size(vid); % Used often

% Get intial orientation
v1 = vid(:,:,1); % initial frame
v2 = medfilt2(imopen(v1,strel('disk',15)),[7 7]);
% v3 = im2bw(v2,graythresh(v2));
v3 = imbinarize(v2);
L = regionprops(v3,'Orientation'); % initial angle

% refangle = L.Orientation;
[yy,~] = find(v3==1);
cent_y = mean(yy);
if cent_y<=(dim(1)/2) % fly pointing up
    refangle = L.Orientation;
else % fly pointing down
    refangle = 180 + L.Orientation;
end

trf     = cell(dim(3),1); % store 2D affine transformations here
fixed   = double(squeeze(vid(:,:,1)));
sz      = imref2d(size(fixed));

% Register each frame with respect to the first frame
tic
for kk = 1:dim(3)
    fprintf([int2str(kk) '\n']);
    z = double(vid(:,:,kk)); % take one frame
    if kk==1
        trf{kk} = imregtform(z,fixed,'rigid',optimizer,metric);
    else
        for jj = kk-1:-1:1
            if(isRigid(trf{jj}))
                break;
            end
        end
        trf{kk} = imregtform(z,fixed,'rigid',optimizer,metric,...
                            'InitialTransformation',trf{jj});
                        %,'PyramidLevels',5);
    end
    reg = imwarp(z,trf{kk},'OutputView',sz);
    vid(:,:,kk) = imrotate(reg,90-refangle,'crop');
    fixed = (fixed*kk + reg)/(kk+1);
end
disp('DONE')
toc
end