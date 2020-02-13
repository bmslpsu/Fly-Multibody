function [hAngles,cPoint,validity,ROI,initframe,finalframe] = headtracker(vid, npoint, playback, showpoint)
% HeadTracker: tracks insect head movments in a rigid terher
%
% Tracks feature on fly head (ususally antenna) & calculates the angle with 
% respect to a specififed center point. Kanade–Lucas–Tomasi feature tracker.
% https://www.mathworks.com/help/vision/ref/vision.pointtracker-system-object.html
%
% Sign convention for angle outputs: [CW = + , CCW = -]
%
% INPUTS:
% 	vidData     :   4D video matrix
%	npoint      :   # of points for tracker
%   playback    :   playback rate (show a frame in increments of "playback")
%                 	If false, then don't show anything (default = 1)
%   showpoint  	:   logical >>> true = debug mode
%
% OUTPUTS:
%   hAngles     :   head angles
%   cPoint      :   head rotation center point
%   validity   	:   validity matrix (nframe x npoint)
%   ROI         :   rectangualr reigon in image used to fid initial
%                   features
%   initframe  	:   first frame of display
%   finalframe 	:   final frame of display
%

if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

vid = squeeze(vid); % remove singleton dimensions
[~,~,nframe] = size(vid); % get size of video

trackFrame = vid(:,:,1); % get 1st frame to start tracker
displayFrame = vid(:,:,1); % get 1st frame for display

% Prompt user to define feature detection area & centerline
fig = figure; clf
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [8 7];
title('Pick area of interest & draw head midline', 'FontWeight', 'bold', 'FontSize',12)
xlabel('Click to continue', 'FontWeight', 'bold', 'FontSize',12)
ax(1) = subplot(1,1,1) ; axis image ; hold on
    imshow(displayFrame) % show frame
    roi = drawrectangle(ax(1)); % draw box around tracking point (antenna)
    centline = drawline(ax(1)); % draw line through head rotation point and head midline
    pause % wait for click to continue
    ROI = round(roi.Position); % get head ROI to detect feature points
    cLine = (centline.Position); % get head center line
close(fig) % close figure

% Coordinates for head rotation point
cPoint.Y = max(cLine(:,2)); % lower y-coordinate is the neck joint
cPoint.X = cLine((cLine(:,2) == cPoint.Y)); % get corrsponding x-coordinate

% Coordinates for initial angle of head midline
mPoint.Y = min(cLine(:,2)); % higher y-coordinate is the midline point
mPoint.X = cLine((cLine(:,2) == mPoint.Y)); % get corrsponding x-coordinate

% Calculate initial angle of head midline
initAngle.head = rad2deg( atan2( mPoint.X - cPoint.X , -(mPoint.Y - cPoint.Y) ) );

% Setup tracking
points	= detectMinEigenFeatures(trackFrame,'ROI',ROI); % detect features
points  = points.selectStrongest(npoint); % get strongest points
tracker = vision.PointTracker('MaxBidirectionalError',5,'NumPyramidLevels',7,...
    'BlockSize',[31 31],'MaxIterations',40); % create tracker object

initialize(tracker,points.Location,trackFrame); % start tracker

% Calculate initial angle to feature (antenna)
[points,validity] = tracker(trackFrame);
pointFrame_disp = insertMarker(displayFrame,points(validity, :),'+');

initFeat.Y = round(mean(points(:,2)));
initFeat.X = round(mean(points(:,1)));

initAngle.feature = rad2deg(atan2(initFeat.X - cPoint.X , -(initFeat.Y - cPoint.Y) ));

% Calculate offset angle from feature to head midline
offsetAngle = initAngle.feature - initAngle.head;

% Show points and lines
fig = figure; clf
ax(1) = subplot(1,1,1) ; axis image ; hold on
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [6 5];
title('Detected Interest Points', 'FontWeight', 'bold', 'FontSize',12)
xlabel('Click to continue', 'FontWeight', 'bold', 'FontSize',12)
    imshow(pointFrame_disp)
    line([cPoint.X ,  mPoint.X ],[cPoint.Y , mPoint.Y],   'Color', 'r', 'Linewidth', 1.5)
    line([cPoint.X , initFeat.X],[cPoint.Y , initFeat.Y], 'Color', 'c', 'Linewidth', 1.5)
pause % wait for click to continue
close(fig) % close figure

% Create display
fig(1) = figure (101); clf
set(fig, 'Visible', 'off')
fColor = 'k'; % figure and main axis color
aColor = 'w'; % axis line color
set(fig,'Color',fColor,'Units','inches','Position',[2 2 9 7])
fig(1).Position(3:4) = [9 7];
% movegui(fig,'center')
figure (101)
    % Raw image with tracking window
    ax(1) = subplot(3,4,[2,3,6,7]); hold on ; cla ; axis image

    % Head angle window
    ax(2) = subplot(3,4,9:12); hold on ; cla ; xlim([0 nframe])
        xlabel('Frame')
        ylabel('Angle (°)')
        h.hAngle = animatedline(ax(2), 'Color', 'b', 'LineWidth', 1);
        % ylim(5*[-1 1])
        
set(ax, 'Color', fColor, 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
    'YColor', aColor, 'XColor',aColor)

% Preallocate vectors to store tracked points & angles
Pos.X    = zeros(nframe,1);
Pos.Y    = zeros(nframe,1);
hAngles  = zeros(nframe,1); 
POINTS   = cell(nframe,1);
validity = zeros(nframe,npoint);

% Track features & calculate head angle for every frame
tic
disp('Tracking...')
for idx = 1:nframe
	% Display frame count every every 100 frames
    if (idx==1) || ~mod(idx,100) || (idx==nframe)
        fprintf('%i\n',idx)
    end
    
    % Get frame & track features
    frame = vid(:,:,idx); % get frame
    [points,validity(idx,:)] = tracker(frame); % detect features
    POINTS{idx} = points; % store points in array
    
	% Calculate average (x,y) position of features
	Pos.Y(idx) = mean(points(:,2));
    Pos.X(idx) = mean(points(:,1));
    
 	% Calculate head midline angle
	hAngles(idx) = rad2deg( atan2( Pos.X(idx) - cPoint.X , -(Pos.Y(idx) - cPoint.Y) ) ) - offsetAngle;
    
    if playback || idx==1 || idx==nframe
        set(fig, 'Visible', 'on')
        % Display image
        if (idx==1) || (~mod(idx,playback)) || (idx==nframe) % at playback rate
            % Show images with tracking annotation
            ax(1) = subplot(3,4,[2,3,6,7]); cla % frame & tracking
            if showpoint % show all tracked points and point mean
                pointFrame = insertMarker(frame, points, '+'); % add points to image
                imshow(pointFrame) % frame with tracked points
                
            	line([cPoint.X , Pos.X(idx) ], [cPoint.Y , Pos.Y(idx)], ... % update line drawn to features
                        'Color', 'c', 'LineWidth', 1)
            else % just show the head midline
                imshow(frame)
            end
                    
            line([cPoint.X ,  Pos.X(idx) - (initFeat.X - mPoint.X)], ... % update line drawn to head midline
                 [cPoint.Y ,  Pos.Y(idx) - (initFeat.Y - mPoint.Y)], ...
                    'Color', 'b', 'LineWidth', 1.5, 'Marker', '.')
        end
    
        % Display angle
        ax(2) = subplot(3,4,9:12); % angles
            addpoints(h.hAngle, idx, hAngles(idx))
            
        if idx==1
            initframe = getframe(fig);
            initframe = initframe.cdata;
        elseif idx==nframe
            finalframe = getframe(fig);
            finalframe = finalframe.cdata;
        end
        
        pause(0.0005) % give time for images to display
    end
end
toc

end