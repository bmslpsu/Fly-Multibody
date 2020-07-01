function [norm_ang,imgstats,initframe] = bodytracker(vid, playback, head_debug)
%% bodytracker: tracks the body angle of an insect in a magnetic tether
%
% Fits an ellipse to the 'on' reigon in each frame (insect body). Blurs the
% image to get rid of head movements and wraps the angle. Can debug by
% displaying a tracking animation.
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
%   INPUT:
%       vid         :   input video matrix
%       playback   	:   playback rate (show a frame in increments of "playback").
%                       If false, then only show the 1st frame. (default = 1)
%       head_debug  :   always bring up the heading angle check-window if true
%                           
%   OUTPUT:
%       norm_ang 	:   normalized & unwrapped body angle [°]
%       imgstats 	:   structure containing some basic image statistics (orientation, centroid, etc.)
%       initframe   :   initial (1st) frame image
%

if nargin < 3
    head_debug = false; % default
    if nargin < 2
        playback = 1; % default
    end
end

if isempty(playback)
    playback = 1; % default
end

if ~rem(playback,1)==0
    warning('Warning: "playback" rounds to nearest integer')
    playback = round(playback);
end

vid = squeeze(vid); % get rid of singleton dimension
vid = flipud(vid); % flip video to arena reference frame
[yp,xp,nframe] = size(vid);  % get size & # of frames of video

% Use the inital frame to find the heading
[~,flip] = findheading(vid(:,:,1), head_debug);
disp('Heading found.')

% Preprocess raw video
disp('Video preprocessing...')
bnvid = false(yp,xp,nframe); % stores the processed video
SE_erode = strel('disk',8,4); % erosion mask
SE_dilate = strel('disk',12,4); % dilation mask
tic
parfor idx = 1:nframe
    % disp(idx)
    frame = vid(:,:,idx); % get raw frame
    bnframe = imbinarize(frame); % binarize
    bnframe = imerode(bnframe, SE_erode); % erode
    bnframe = imdilate(bnframe, SE_dilate); % dilate
    bnvid(:,:,idx) = logical(bnframe); % store bianary frame
end
toc

pause(1)
close all

% Set some parameters
offset  = 90;  % shift the reference frame so 0° is the top vertical axis in video [°]
dthresh = 160; % threshold for detecting >180° flips in ellipse orientation, or angles > 360°
shift   = 0;   % shift to keep angle wrapped [°] (dynamic)

% Create display
fig(1) = figure (100); clf
set(fig, 'Visible', 'off')
fColor = 'k'; % figure and main axis color
aColor = 'w'; % axis line color
set(fig, 'Color', fColor, 'Units', 'inches', 'Position', [2 2 9 7])
fig(1).Position(3:4) = [9 7];
% movegui(fig,'center')
figure (100)
    % Raw image window
    ax(1) = subplot(3,2,[1,3]); hold on ; cla ; axis image

    % Processed video window
    ax(2) = subplot(3,2,[2,4]); hold on ; cla ; axis image

    % Body angle window
    ax(3) = subplot(3,2,5:6); hold on ; cla ; xlim([0 nframe])
        xlabel('Frame')
        ylabel('Angle (°)')
        h.raw_angle  = animatedline(ax(3),'Color','b','LineWidth',1); % debugging
        h.norm_angle = animatedline(ax(3),'Color','r','LineWidth',1);

set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
    'YColor', aColor, 'XColor',aColor)

% Preallocate vectors to store tracked angles
raw_ang  = nan(nframe,1); % stores raw angles calulated by ellipse fit [°]
norm_ang = nan(nframe,1); % stores normalized/unwrapped angles [°]

tic
disp('Tracking...')
for idx = 1:nframe
    % Display frame count every every 100 frames
    if (idx==1) || ~mod(idx,100) || (idx==nframe)
        fprintf('%i\n',idx)
    end
    
    % Get images
    frame = vid(:,:,idx); % raw frame
    bnframe = bnvid(:,:,idx); % processed frame
    
  	% Calculate angle
    tempstats = regionprops(bnframe,'Centroid','Area','BoundingBox','Orientation', ...
        'MajorAxisLength','MinorAxisLength'); % image reigon properties
    
    imgstats(idx) = tempstats(1); % get stats only for largest feature
    
    centroid = imgstats(idx).Centroid; % centroid of image
    L = imgstats(idx).MajorAxisLength / 2; % long axis of image
    
    raw_ang(idx)  = (imgstats(idx).Orientation - offset); % raw angle [°]
    norm_ang(idx) = -(imgstats(idx).Orientation - offset + shift); % normalized, unwrapped angle [°]
 	
    % Check for changes in angle > 180°. Correct for the ellipse fit and unwrap angles.
    if idx > 1
        dbody = norm_ang(idx) - norm_ang(idx-1); % change in body angle between frames [°]
        magd = abs(dbody); % magnitude of change [°]
        signd = sign(dbody); % direction of change
        if magd > dthresh % 180 or more flip
            flipn = round(magd/180); % how many 180° shifts we need
         	shift = -signd*flipn*180; % shift amount [°]
            norm_ang(idx) = norm_ang(idx) + shift; % normalized, unwrapped angle [°]
        end
    elseif idx==1
        % Flip heading by 180° if the heading is in the wrong direction
        if flip
            norm_ang(idx) = norm_ang(idx) + 180;
        end
        
       	% Make start angle positive
        if norm_ang(idx) < 0
            norm_ang(idx) = norm_ang(idx) + 360;
        end
    end

    if playback || idx==1
        set(fig, 'Visible', 'on')
        % Display images
        if (idx==1) || (~mod(idx,playback)) || (idx==nframe) % display at playback rate
            % Get approximate location of head
            head = centroid + L*[sind(norm_ang(idx)), -cosd(norm_ang(idx))];
            heading = [centroid ; head];
            abdomen = centroid - L*[sind(norm_ang(idx)), -cosd(norm_ang(idx))];
            reverse = [centroid ; abdomen];

            % Show images with tracking annotation
            ax(1) = subplot(3,2,[1,3]); cla % raw
                imshow(frame) % frame
                
                h.heading = semi_ellipse(centroid, L, 0.5, 0.90, 180 - norm_ang(idx), 'r');
                delete([h.heading{2:4}])
                alpha(h.heading{1},0.3)
                h.heading{1}.LineStyle = 'none';
                
            ax(2) = subplot(3,2,[2,4]); cla % processed
                imshow(bnframe) % frame

                % Show bounding box
                h.rect = rectangle('Position', imgstats(idx).BoundingBox, 'EdgeColor', 'g', ...
                                                                'LineWidth', 1);

                % Show ellipse
                h.ellps = ellipse(centroid, 2*L, 0.5, 0.90, 180 - norm_ang(idx), 'r');
                delete([h.ellps{[1,3:6]}])
                h.ellps{2}.Color = 'r';
                h.ellps{2}.LineWidth = 1;
                
                h.heading = semi_ellipse(centroid, L, 0.5, 0.90, 180 - norm_ang(idx), 'r');
                delete([h.heading{2:4}])
                alpha(h.heading{1},0.3)
                h.heading{1}.LineStyle = 'none';

                % Show centroid & heading
                plot(heading(:,1), heading(:,2), '-r', 'LineWidth', 1) % centroid
                plot(reverse(:,1), reverse(:,2), '-b', 'LineWidth', 1) % heading
                plot(centroid(1),  centroid(2),  '.g', 'MarkerSize', 20) % reverse heading
        end
    
        % Display angle
        ax(3) = subplot(3,2,5:6);
            % addpoints(h.raw_angle, idx, raw_angle(idx)) % debugging
            addpoints(h.norm_angle, idx, norm_ang(idx))
            
        if idx==1 % get 1st frame
            initframe = getframe(fig);
            initframe = initframe.cdata;
        end
        
        pause(0.0005) % give time for images to display
    end
end
toc
end