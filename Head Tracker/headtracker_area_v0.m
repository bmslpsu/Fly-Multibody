function [hAngles,cPoint,validity,ROI,initframe,finalframe] = headtracker_area(vid, playback, showpoint)
%% headtracker_area: tracks insect head movments in a rigid terher
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

playback = 1;

if ~rem(playback,1)==0
    warning('Warning: "playback" roundes to nearest integer')
    playback = round(playback);
end

vid = squeeze(vidData); % remove singleton dimensions
[yp,xp,n_frame] = size(vid); % get size of video

% [bAngle] = bodytracker(frame, false, true, false);

%% Get image stats across all frames
clc
filt_vid = vid;
clear all_stats
for n = 1:n_frame
    A = cell(1);
    A{1} = vid(:,:,n);
    A{end+1} = medfilt2(A{end} ,[3 3]);
	A{end+1} = (imbinarize(A{end}));
    A{end+1} = uint8(bwareaopen(A{end},50));
    filt_vid(:,:,n) = A{end};
    
    stats = regionprops(A{end},'Centroid','MajorAxisLength','MinorAxisLength','BoundingBox');
    all_stats(n) = stats;
end

%% Find bounding box across all frames
bounding_box = median(cat(1,all_stats.BoundingBox),1);
flex = [0 0.1];
cut = 0.3;
x1 = floor((1 - flex(1)) * bounding_box(1));
y1 = floor((1 - flex(2)) * bounding_box(2));
x2 = ceil((1 + flex(1)) * (bounding_box(1) + bounding_box(3)));
y2 = ceil((1 - cut) * (bounding_box(2) + bounding_box(4)));
xr = x1:x2;
yr = y1:y2;
out_range_x = (xr > xp) | (xr < 1);
out_range_y = (yr > yp) | (yr < 1);
xr = xr(~out_range_x);
yr = yr(~out_range_y);

%% Cut, filter, and morph image to find outline
se = strel('disk',2);
cut_vid = vid(yr,xr,:);
cut_filt_vid = false(size(cut_vid));
clear cut_stats
for n = 1:n_frame
    A = cell(1);
    A{1} = cut_vid(:,:,n);
    A{end+1} = medfilt2(A{end}, [5 5]);
	A{end+1} = imbinarize(A{end});
    A{end}(end,:) = true;
    A{end+1} = imfill(A{end}, 'holes'); 
    A{end+1} = bwareaopen(A{end},50);
    
    cent_frame = imdilate(A{end}, se);
    stats = regionprops(cent_frame,'Centroid');
    cut_stats(n) = stats;
    
    A{end+1} = imerode(A{end}, se);
    A{end+1} = imfill(A{end}, 'holes'); 
    A{end+1} = bwmorph(A{end}, 'bridge',3);
    A{end+1} = bwmorph(A{end},'remove');
    %montage(A)

    cut_filt_vid(:,:,n) = logical(A{end});
 	%pause()
end
cent = median(cat(1,cut_stats.Centroid),1);
dim = size(cut_filt_vid);

%%
clear ax
roll_cal = 36.33;
playback = 0;
export = false;
if export
    Fs = 50;
    targetdir = 'H:\EXPERIMENTS\RIGID\Experiment_SOS_v2\movie';
    fname = 'Test_SOS_v6.mp4';
    VID = VideoWriter(fullfile(targetdir,fname),'MPEG-4');
    VID.FrameRate = Fs;
    open(VID)
end

if playback
    fig = figure (10) ; clf
    set(fig, 'Color', [0 0 0.1], 'Units', 'inches', 'Position', [4 1 10 9])
        ax(1) = subplot(4,3,2); cla ; hold on ; axis image
        ax(2) = subplot(4,3,5); cla ; hold on ; axis image
        ax(3) = subplot(4,3,1); cla ; hold on ; axis image
        ax(4) = subplot(4,3,3); cla ; hold on ; axis image
        ax(5) = subplot(4,3,4); cla ; hold on ; axis image
        ax(6) = subplot(4,3,6); cla ; hold on ; axis image
        ax(7) = subplot(4,3,7); cla ; hold on ; axis image
        ax(8) = subplot(4,3,8); cla ; hold on ; axis image
        ax(9) = subplot(4,3,9); cla ; hold on ; axis tight ; ylabel('Intensity')
        ax(10) = subplot(4,3,10:12); cla ; hold on ; ylabel('Angle (°)') ; xlabel('Frame')
            h.yaw = animatedline(ax(10), 'Color', 'm', 'LineWidth', 1);
            h.roll = animatedline(ax(10), 'Color', 'r', 'LineWidth', 1);
            xlim([1 n_frame])

        set(ax(1:10), 'XColor', 'w', 'YColor', 'w', 'Box', 'on', 'FontSize', 9)
        set(ax(1:8), 'Color', 'none')
        set(ax(9), 'Color', 'k')
        set(ax(10), 'Color', 'k', 'LineWidth', 1, 'Box', 'off')
        set(ax(1:9), 'XTick', [1 20:20:dim(2)], 'XLim', [0 dim(2)])
        set(ax(1:8), 'YTick', [1 20:20:dim(1)]) 
end

hAngles = nan(n_frame,4);
tic
for n = 1:n_frame
    disp(n)
    raw = cut_vid(:,:,n);
    frame = cut_filt_vid(:,:,n);
    
    left = frame;
    left(:,round(cent(1)):dim(2)) = false;
    right = frame;
    right(:,1:round(cent(1))) = false;

    xx = zeros(dim(1),2);
    xR = zeros(dim(1),1);
    xL = zeros(dim(1),1);
    for r = 1:dim(1)
        all_x = find(frame(r,:));
        if length(all_x) > 1
            xx(r,:) = [all_x(1) all_x(end)];
        end
        x_right = find(right(r,:), 1, 'last');
        x_left = find(left(r,:), 1, 'first');
        if ~isempty(x_right)
            xR(r) = find(right(r,:), 1, 'last');
        end
        if ~isempty(x_left)
            xL(r) = find(fliplr(left(r,:)), 1, 'last');
        end
    end
    [b, a] = butter(3, 0.5, 'low');
    xR_filt = filtfilt(b, a, xR);
    xL_filt = filtfilt(b, a, xL);

    [neck_right] = calculateHeadAngle(xR_filt, right, false);
    [neck_left] = calculateHeadAngle(xL_filt, left, false);
    
    neck_all.x = [neck_left.x ; neck_right.x];
    neck_all.y = [neck_left.y ; neck_right.y];
    
    neck_all.x = neck_all.x(~isnan(neck_all.x));
    neck_all.y = neck_all.y(~isnan(neck_all.y));
    
  	c = polyfit(neck_all.x, neck_all.y, 1);
    neck_all.x_test = (1:dim(2))';
    neck_all.y_test = polyval(c, neck_all.x_test);
    neck_all.angle = rad2deg(atan(c(1)));
    
    cent_offset = 0.2*min([neck_left.width, neck_right.width]);
    cPoint.x = cent(1);
    cPoint.y = nanmean([neck_left.peak, neck_right.peak]) - cent_offset;
    
    [~,intersect.x] = min(abs(neck_all.x_test - cent(1)));
    intersect.y = neck_all.y_test(intersect.x);

    hAngles(n,1) = neck_right.angle;
    hAngles(n,2) = neck_left.angle;
    hAngles(n,3) = neck_all.angle;
    yaw = hAngles(n,3);
    
    rot_frame = imrotate(raw, hAngles(n,3), 'bilinear', 'crop');
    rot_frame = medfilt2(rot_frame,5*[1 1]);
    rot_frame_bw = imbinarize(rot_frame);
    all_on = find(max(rot_frame_bw,[],2), 1, 'first');
    %rot_frame_cut = rot_frame;
    %rot_frame_cut(round(cPoint.y):dim(1),:) = 0;
    rot_frame_head = rot_frame(all_on:round(intersect.y),:);
    
    test(n) = size(rot_frame_head,1);
    
    stats = regionprops(imbinarize(rot_frame_head),'Centroid','BoundingBox');
    midI = round(stats(1).Centroid(2));
    span = round(0.35*stats(1).BoundingBox(4) / 2);
    eye_range = midI-span:midI+span;
    eye_box.x = [1 1 dim(1) dim(1)];
    eye_box.y = [eye_range(1) eye_range(end) eye_range(end) eye_range(1)];
    eye_intensity = double(rot_frame_head(eye_range,:));
    eye_intensity_mean = median(eye_intensity,1);
    
    [pks,locs,w,p] = findpeaks(eye_intensity_mean, 'MinPeakHeight', 100, 'MinPeakDistance', 5, ...
                                    'MinPeakWidth', 5, 'MinPeakProminence', 5);
    eye_peaks = locs([1,end]);
    eye_peak_int = pks([1,end]);
    peak_int_mean = min(0.95*eye_peak_int);
    
    [xi,~] = polyxpoly([1 dim(2)],[peak_int_mean peak_int_mean],1:dim(2),eye_intensity_mean);
    eye_left = [xi(1) , xi(2)];
    eye_right = [xi(end-1) , xi(end)];
    eye_left_width = diff(eye_left);
    eye_right_width = diff(eye_right);
    roll_idx = (eye_left_width - eye_right_width) / (eye_left_width + eye_right_width);
    roll = roll_cal * roll_idx;
    
    hAngles(n,4) = roll;
                            
    r = stats(1).BoundingBox(4);
    rpoint.x = intersect.x + r*sind(hAngles(n,3));
    rpoint.y = intersect.y - r*cosd(hAngles(n,3));
 	
    if playback && ( (n==1) || (~mod(n,playback)) || (n==n_frame) )
        ax(1) = subplot(4,3,2); cla ; hold on
            imshow(raw)
            plot([cent(1) cent(1)], [1 dim(1)], '--c', 'LineWidth', 0.5);
            plot([1 dim(2)], [intersect.x intersect.y], '--c', 'LineWidth', 0.5);
            plot([intersect.x rpoint.x], [intersect.y rpoint.y], 'm', 'LineWidth', 2)
            title(['Yaw = ' num2str(yaw) '°'], 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 14)
            
            plot(neck_all.x, neck_all.y, '.b', 'MarkerSize', 5)
            plot(neck_all.x_test, neck_all.y_test, 'm', 'LineWidth', 2)
            plot(cPoint.x, cPoint.y, '.g', 'MarkerSize', 15)
            plot(intersect.x, intersect.y, '.m', 'MarkerSize', 15)

        ax(2) = subplot(4,3,5); cla ; hold on
            imshow(255*frame)
            plot([cent(1) cent(1)], [1 dim(1)], '--c', 'LineWidth', 0.5);
            plot([1 dim(2)], [intersect.x intersect.y], '--c', 'LineWidth', 0.5);
            plot([intersect.x rpoint.x], [intersect.y rpoint.y], 'm', 'LineWidth', 1)
            
            plot(neck_all.x, neck_all.y, '.b', 'MarkerSize', 5)
            plot(neck_all.x_test, neck_all.y_test, 'm', 'LineWidth', 1)
            plot(cPoint.x, cPoint.y, '.g', 'MarkerSize', 15)
            plot(intersect.x, intersect.y, '.m', 'MarkerSize', 15)

        ax(3) = subplot(4,3,1); cla ; hold on
            imshow(0*frame)
            plot(flipud(xL_filt), dim(1):-1:1, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
            plot(neck_left.dx, neck_left.peak, '.g', 'MarkerSize', 15)

        ax(4) = subplot(4,3,3); cla ; hold on
            imshow(0*frame)
            plot(flipud(xR_filt), dim(1):-1:1, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
            plot(neck_right.dx, neck_right.peak, '.g', 'MarkerSize', 15)

        ax(5) = subplot(4,3,4); cla ; hold on
            imshow(left)
            patch([1 1 dim(2) dim(2)], [neck_left.edge_range(1) neck_left.edge_range(end) ...
                neck_left.edge_range(end) neck_left.edge_range(1)], 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
            plot(xL,1:dim(1),'m', 'LineWidth', 1)
            plot(neck_left.dx, neck_left.peak, '.g', 'MarkerSize', 15)
            plot(neck_left.x, neck_left.y, '.b', 'MarkerSize', 5)
            plot(neck_left.x_test, neck_left.y_test, 'y', 'LineWidth', 1)

        ax(6) = subplot(4,3,6); cla ; hold on
            imshow(right)
            patch([1 1 dim(2) dim(2)], [neck_right.edge_range(1) neck_right.edge_range(end) ...
                neck_right.edge_range(end) neck_right.edge_range(1)], 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
            plot(xL,1:dim(1),'m', 'LineWidth', 1)
            plot(neck_right.dx, neck_right.peak, '.g', 'MarkerSize', 15)
            plot(neck_right.x, neck_right.y, '.b', 'MarkerSize', 5)
            plot(neck_right.x_test, neck_right.y_test, 'y', 'LineWidth', 1)
            
        ax(7) = subplot(4,3,7); cla ; hold on
            imshow(rot_frame)
            plot([1 dim(2)], [cPoint.y cPoint.y], '--g', 'LineWidth', 0.5)
            plot([1 dim(2)], [all_on all_on], '--g', 'LineWidth', 0.5)
            plot(cPoint.x, cPoint.y, '.g', 'MarkerSize', 15)
            
        ax(8) = subplot(4,3,8); cla ; hold on
            imshow(rot_frame_head)
            patch(eye_box.x, eye_box.y, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
          	plot(eye_left, [midI midI],  '-r', 'LineWidth', 3)
            plot(eye_right, [midI midI],  '-r', 'LineWidth', 3)
            title(['Roll = ' num2str(roll) '°'], 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 14)
            axis on
            
        ax(9) = subplot(4,3,9); cla ; hold on
            plot(eye_intensity', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            plot(eye_intensity_mean, 'w', 'LineWidth', 2)
            plot(eye_peaks, eye_peak_int, '.g', 'MarkerSize', 20)
            plot([1 dim(2)], [peak_int_mean peak_int_mean], '-b', 'LineWidth', 1)
            plot(eye_left, [peak_int_mean peak_int_mean],  '-r', 'LineWidth', 3)
            plot(eye_right, [peak_int_mean peak_int_mean],  '-r', 'LineWidth', 3)
            
            ylim([100 260])
            ax(9).Position([1,3]) = ax(4).Position([1,3]);
    end
    
    if playback
     	ax(10) = subplot(4,3,10:12); hold on
            addpoints(h.yaw, n, yaw)
            addpoints(h.roll, n, roll)
        set(ax, 'Visible', 'on')
    end
    
    if export && playback
     	fig_frame = getframe(fig);
    	writeVideo(VID,fig_frame);
    end
    pause(0.001)
end
toc
if export && playback
    close(VID)
end
end

function [neck] = calculateHeadAngle(edge, I, debug)
%% calculateHeadAngle:
%   INPUTS:
%       edge    : far edge for each row
%       I       : left or right image
%       debug   : show plots
%
%   OUTPUTS:
%       neck    : neck structure
%
dim = size(I);
[pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                                'MinPeakProminence', 10, 'SortStr', 'descend');
idx = (1:length(edge))';
if isempty(locs)
    warning('Small peak')
    [pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                           	'MinPeakProminence', 5);
    if isempty(locs)
        warning('Can''t find angle')
        pause(0.1)
        neck.peak = nan;
        neck.dx = nan;
        neck.width = nan;
        neck.edge_range = nan;
        neck.edge = nan;
        neck.image = nan;
        neck.image_all = nan;
        neck.angle = nan;
        neck.x = nan;
        neck.y = nan;
        neck.x_test = nan;
        neck.y_test = nan;
        return
    else
        locs = locs(end);
        pks = pks(end);
    end
end

neck.peak = locs(1);
neck.dx = -pks(1);
neck.width = w(1);

offset = 0.2*neck.width;

neck.edge_range = ((neck.peak  - ceil(neck.width/1.5) ):floor(neck.peak - offset))';
neck.edge = (edge(neck.edge_range));
neck.image = I(neck.edge_range,:);
neck.image_all = false(size(I));
neck.image_all(neck.edge_range,:) = neck.image;

[y, x] = find(neck.image_all);
neck.x = x;
neck.y = y;
neck.x_test = 1:size(I,2);

c = polyfit(x, y, 1);

neck.y_test = polyval(c, neck.x_test);
hangle = rad2deg(atan(c(1)));
neck.angle = hangle;

if debug
    fig = figure (109);
    set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 4 6])
    fig.Position
    ax(1) = subplot(1,1,1); hold on; axis image
        imshow(255*I, 'InitialMagnification', 400) ; hold on
        patch([1 1 dim(2) dim(2)], [neck.edge_range(1) neck.edge_range(end) ...
            neck.edge_range(end) neck.edge_range(1)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        plot(edge,idx,'m', 'LineWidth', 1)
        plot(neck.dx, neck.peak, '.g', 'MarkerSize', 15)
        plot(neck.x, neck.y, '.b', 'LineWidth', 1)
        plot(neck.x_test, neck.y_test, '-r', 'LineWidth', 1)
    set(ax, 'Box', 'on')
end

end