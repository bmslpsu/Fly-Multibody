function [vid,trf] = regvid(vid)
%% register_video: register each frame of a video to the same reference
%
%   INPUT:
%       vid   	:   raw video
%
%   OUTPUT:
%       vid   	:   registered video
%       trf   	:   2D affine transformation for each frame
%
clearvars -except vidData t_v

vid = vidData;
vid = squeeze(vid);

%
% Set up optimizer
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 1e-3;
optimizer.MaximumIterations = 150;
optimizer.MaximumStepLength = 0.03;
optimizer.MinimumStepLength = 0.0002;

% vid = flipud(squeeze(vid)); % In case it's 4D
dim = size(vid); % Used often

% Get intial orientation
v1 = vid(:,:,1); % initial frame
[refangle,flip,fixed] = findheading(v1, 1, 1/4);
if flip
    refangle = refangle - 180;
end

trf = cell(dim(3),1); % store 2D affine transformations here
sz = imref2d(size(fixed));

% z = vid(:,:,900);
% moving_reg = imregister(z, fixed, 'rigid', optimizer, metric,'InitialTransformation', init_trf);
% montage({fixed, z, flipud(moving_reg)})

%% Register each frame with respect to the first frame
tic
trf{1} = imregtform(vid(:,:,1), fixed, 'rigid', optimizer, metric);
regvid = uint8(zeros(dim));
regvid(:,:,1) = imwarp(vid(:,:,1), trf{1}, 'OutputView', sz);
for kk = 2:dim(3)
    fprintf([int2str(kk) '\n']);
    z = vid(:,:,kk); % take one frame
    trf{kk} = imregtform(z, fixed, 'rigid', optimizer, metric, 'InitialTransformation', trf{kk-1});                    
    reg = imwarp(z, trf{kk}, 'OutputView', sz);
    %regvid(:,:,kk) = imrotate(reg, refangle,'crop');
    regvid(:,:,kk) = reg;
    
    %fixed = (fixed*kk + reg)/(kk+1);
    %cla
    %montage({fixed, z, reg,})
    %pause
end
disp('DONE')
toc
end