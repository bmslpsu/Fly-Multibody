function [norm_ang,imgstats] = bodytracker(vid, playback)
%% bodytracker: tracks the body angle of an insect in a magnetic tether
%
%   INPUT:
%       vid         :   input video matrix
%       playback   	:   playback rate (show a frame in increments of "playback").
%                       If false, then don't show anything.
%                           default = 1
%   OUTPUT:
%       norm_ang 	:   normalized, unwrapped angle [°]
%       imgstats 	:   structure containing some basic image statistics (orientation, centroid, etc.)
%

if nargin < 2
    playback = 1; % default
end

vid = flipvid(vid,'lr'); % flip video to arena reference frame
[yp,xp,nframe] = size(vid);  % get size & # of frames of video

% Preprocess raw video
disp('Video preprocessing...')
bnvid = false(yp,xp,nframe); % stores the processed video
SE_erode = strel('disk',8,4);   % erosion mask
SE_dilate = strel('disk',12,4); % dilation mask
tic
parfor idx = 1:nframe
    % disp(idx)
    frame = vid(:,:,idx);
    bnframe = imbinarize(frame);
    bnframe = imerode(bnframe,  SE_erode);
    bnframe = imdilate(bnframe, SE_dilate);
    bnvid(:,:,idx) = logical(bnframe);
end
toc

% Use 1st image to find intial heading
first_stats = regionprops(bnvid(:,:,1000),'Centroid','Area','BoundingBox','Orientation','Image', ...
    'MajorAxisLength','MinorAxisLength'); % image reigon properties

init_frame = first_stats.Image;
check_frame = imrotate(init_frame, 90 - first_stats.Orientation, 'loose');
check_frame = imerode(check_frame, strel('disk',20,8));

top_ratio = sum(check_frame(1:fix(numel(check_frame)/2))) ./ ...
                sum(check_frame(fix(numel(check_frame)/2)+1:numel(check_frame)));
            
if top_ratio < 1
    flip = true
else
    flip = false
end
            
% figure (1) ; clf
% montage({init_frame,check_frame})

% Create display
if playback
    fig(1) = figure; clf
    fColor = 'k'; % figure and main axis color
    aColor = 'w'; % axis line color
    set(fig,'Color',fColor,'Units','inches','Position',[2 2 9 7])
    fig(1).Position(3:4) = [9 7];
    % movegui(fig,'center')
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
end

raw_ang  = nan(nframe,1); % stores raw angles calulated by ellipse fit [°]
norm_ang = nan(nframe,1); % stores normalized/unwrapped angles [°]

offset = 90; % shift the reference frame so 0° is the top vertical axis in video [°]
dthresh = 160; % threshold for detecting 180 flips in ellipse orientation, or angles > 360 [°]
shift = 0; % shift to keep angle wrapped [°] (dynamic)

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
    imgstats(idx) = regionprops(bnframe,'Centroid','Area','BoundingBox','Orientation','Image', ...
        'MajorAxisLength','MinorAxisLength'); % image reigon properties
    
    centroid = imgstats(idx).Centroid;
    L = imgstats(idx).MajorAxisLength / 2;
    
    raw_ang(idx)  = -(imgstats(idx).Orientation + offset); % raw angle [°]
    norm_ang(idx) = -(imgstats(idx).Orientation + offset + shift); % normalized, unwrapped angle [°]
    
    if flip
        norm_ang(idx) = norm_ang(idx) + 180;
    end
 	
    % Check for changes in angle > 180 degree. correct for the ellipse
    % fit and unwrap angles
    if idx > 1
        dbody = norm_ang(idx) - norm_ang(idx-1); % change in body angle between frames [°]
        magd = abs(dbody); % magnitude of change [°]
        signd = sign(dbody); % direction of change
        if magd > dthresh % 180 or more flip
            flipn = round(magd/180); % how many 180 shifts we need
         	shift = -signd*flipn*180; % shift amount [°]
            norm_ang(idx) = norm_ang(idx)  + shift; % normalized, unwrapped angle [°]
        end
    elseif idx==1 % start angle as positive
        if norm_ang(idx) < 0
            norm_ang(idx) = norm_ang(idx) + 360;
        end
    end

    if playback
        % Display images
        if (idx==1) || (~mod(idx,playback)) || (idx==nframe) % display at playback rate
            % Get approximate location of head
            head = centroid + L*[sind(norm_ang(idx)), -cosd(norm_ang(idx))];
            heading = [centroid ; head];
            abdomen = centroid - L*[sind(norm_ang(idx)), -cosd(norm_ang(idx))];
            reverse = [centroid ; abdomen];

            % Show images with tracking annotation
            subplot(3,2,[1,3]); cla % raw
                imshow(1.5*frame)
            subplot(3,2,[2,4]); cla % processed
                imshow(bnframe) % frame

                % Show bounding box
                h.rect = rectangle('Position', imgstats(idx).BoundingBox, 'EdgeColor', 'g', ...
                    'LineWidth', 1); % bounding box

                % Show ellipse
                h.ellps = ellipse(centroid, 2*L, 0.5, 0.90, -norm_ang(idx), 'r');
                delete([h.ellps{[1,3:6]}])
                h.ellps{2}.Color = 'r';
                h.ellps{2}.LineWidth = 1;
                hold on

                % Show centroid & heading
                plot(heading(:,1), heading(:,2), '-r', 'LineWIdth', 1) % centroid
                plot(reverse(:,1), reverse(:,2), '-b', 'LineWIdth', 1) % heading
                plot(centroid(1),  centroid(2),  '.g', 'MarkerSize', 20) % reverse heading
        end
    
        % Display angle
        subplot(3,2,5:6)
            % addpoints(h.raw_angle,  idx, raw_angle(idx)) % debugging
            addpoints(h.norm_angle, idx, norm_ang(idx))

        pause(0.0005) % give time for images to display
    end
end
toc

end
