classdef headtracker_area_v1
    %% headtracker_area: track head yaw and roll
    %   
    
    properties %(SetAccess = private)
        vid                 % input video data
        dim                 % input video dimensions
        bw_vid              % binarized video
        vid_props           % binarized video properties
        bound_box           % median bounding box
        body_yaw            % body yaw
        
        head_vid            % head frame of input video
        head_bw_vid         % cut bw head video
        head_out_vid        % cut outline of bw video
        head_dim            % cut video dimensions
        head_cent           % median centroid of cut video
        head_stable_vid     % yaw stabilized video in head reference frame
        head_start          % start of head pixels in each frame
        head_start_min      % head start minimum across frames
        head_iso_vid        % isolated head in yaw stabilzed head video
        head_iso_props      % isoltaed video properties
        head_radius         % head radius [pixels]
        head_tip            % head tip point

        neck_left           % neck properties on left side
        neck_right          % neck properties on right side
        neck                % neck properties
        pivot               % pivot point from yaw line intersection with midline
        yaw                 % head yaw
        saturation          % saturation [L, R] (boolean)
        
        eye                 % eye properties
        roll_idx            % head roll ratio
        roll_cal            % head roll calibration factor
        roll                % head roll
    end
    
    methods
        function obj = headtracker_area_v1(vid)
            % headtracker_area: Construct an instance of this class
            %   
            obj.vid = squeeze(vid); % input video
            obj.dim = size(obj.vid); % video size
            
            [obj.vid_props, obj.bw_vid, obj.bound_box, obj.body_yaw] = get_vid_props(obj, vid, [], false);
            [obj] = get_head_vid(obj, [0 0.1], 0.3);
            [~, ~, obj] = get_yaw(obj, 0.5);
            [~, obj] = get_roll(obj, 0.35, 36.33);
        end
        
        function [vid_props, bw_vid, bound_box, body_yaw] = get_vid_props(obj, vid, fill_ratio, debug)
            % get_image_statsL get image stats from input video
            %   
            
            fprintf('Getting video properties ...')
            se = strel('disk',4);
            bw_vid = false(size(vid));
            for n = 1:obj.dim(3)
                %disp(n)
                A = cell(1);
                A{1} = vid(:,:,n);
                A{end+1} = medfilt2(A{end} ,[3 3]);
                A{end+1} = imbinarize(A{end});
                A{end+1} = bwareaopen(A{end},50);
                %A{end+1} = imdilate( imdilate(A{end}, se), se);
                A{end+1} = imdilate(A{end}, se);
                
                image_props = regionprops(A{end},'Centroid','MajorAxisLength','MinorAxisLength',...
                                            'BoundingBox','Orientation','Image');
             	
                if length(image_props) > 1
                    mjaxes = [image_props.MajorAxisLength];
                    [~,idx] = sort(mjaxes, 'descend');
                    idx = idx(1);
                else
                    idx = 1;
                end                  
                vid_props(n) = image_props(idx);
                
                if ~isempty(fill_ratio)
                    bb = vid_props(n).BoundingBox;
                    x1 = floor( bb(1) );
                    y1 = floor( bb(2) + (1-fill_ratio)*bb(4) );
                    x2 = floor( bb(1) + bb(3) );
                    y2 = floor( bb(2) + bb(4) );
                    xr = x1:x2;
                    yr = y1:y2;
                    temp = A{end};
                    temp(yr,xr) = true;
                    A{end+1} = temp;
                    A{end+1} = imfill(A{end}, 'holes');
                    A{end+1} = bwmorph(A{end},'remove');
                end
                
                bw_vid(:,:,n) = A{end};

                if debug
                   montage(A)
                   pause
                end
            end
            bound_box = median(cat(1,vid_props.BoundingBox),1);
            body_yaw = [vid_props.Orientation];
            body_yaw(body_yaw < -40) = body_yaw(body_yaw < -40) + 180;
            body_yaw = body_yaw - 90;
            fprintf(' Done. \n')
        end
        
        function [obj] = get_head_vid(obj, flex, cut)
            % get_head_vid: 
            %   flex: 1x2 proportion to extend frame by in each direction
            %   cut: proportion to cut frame from bottom
            %
            fprintf('Getting head reigon ...')
            bb = obj.bound_box;
            x1 = floor((1 - flex(1)) * bb(1));
            y1 = floor((1 - flex(2)) * bb(2));
            x2 = ceil((1 + flex(1)) * (bb(1) + bb(3)));
            y2 = ceil((1 - cut) * (bb(2) + bb(4)));
            xr = x1:x2;
            yr = y1:y2;
            out_range_x = (xr > obj.dim(2)) | (xr < 1);
            out_range_y = (yr > obj.dim(1)) | (yr < 1);
            xr = xr(~out_range_x);
            yr = yr(~out_range_y);

            obj.head_vid = obj.vid(yr,xr,:);
            obj.head_out_vid = false(size(obj.head_vid));
            se_dilate = strel('disk',8);
            se_erode = strel('disk',2);
            %se_line = strel('line',5,90);
            for n = 1:obj.dim(3)
                A = cell(1);
                A{1} = 2*obj.head_vid(:,:,n);
                A{end+1} = medfilt2(A{end}, [5 5]);
                A{end+1} = imbinarize(A{end});
                A{end}(end,:) = true;
                A{end+1} = imfill(A{end}, 'holes'); 
                A{end+1} = bwareaopen(A{end},300);

                cent_frame = imdilate( imdilate(A{end}, se_dilate), se_dilate);
                image_props = regionprops(cent_frame,'Centroid','MajorAxisLength');               
                if length(image_props) > 1
                    mjaxes = [image_props.MajorAxisLength];
                    [~,idx] = sort(mjaxes, 'descend');
                    idx = idx(1);
                else
                    idx = 1;
                end                  
                cut_props(n) = image_props(idx);
                
                % Get outline of head/neck
                A{end+1} = imerode(A{end}, se_erode);
                A{end+1} = imfill(A{end}, 'holes');
                A{end+1} = bwareaopen(A{end},300);
                %A{end+1} = bwmorph(A{end}, 'bridge',3);
                
                
                obj.head_bw_vid(:,:,n) = A{end};
                
                A{end+1} = bwmorph(A{end},'remove');
                
                obj.head_out_vid(:,:,n) = A{end};
            end
            obj.head_dim = size(obj.head_out_vid);
            obj.head_cent = median(cat(1,cut_props.Centroid),1);
            fprintf(' Done. \n')
        end
        
        function [yaw, pivot, obj] = get_yaw(obj, Fc_n)
            % get_yaw: find yaw angle of head in each frame
            %  Fc_n: normalized filtering coefficent for finding neck joint
            %
            fprintf('Calculating yaw angles ...')
            obj.head_stable_vid = uint8(zeros(size(obj.head_vid)));
            for n = 1:obj.head_dim(3)
                %fprintf('%i \n', n)
                raw = obj.head_vid(:,:,n); % raw for stabilization
                frame = obj.head_out_vid(:,:,n); % frame for analysis
                
                % Get left and right segments of frame
                left = frame;
                left(:,round(obj.head_cent(1)):obj.head_dim(2)) = false;
                right = frame;
                right(:,1:round(obj.head_cent(1))) = false;
                
                % Get the farthest 'on' pixel from center line for left and right images
                xx = zeros(obj.head_dim(1),2);
                xR = zeros(obj.head_dim(1),1);
                xL = zeros(obj.head_dim(1),1);
                for r = 1:obj.head_dim(1)
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
                
                % Filter curves
                [b, a] = butter(3, Fc_n, 'low');
                xR_filt = filtfilt(b, a, xR);
                xL_filt = filtfilt(b, a, xL);
                
                % Get bottom edge of neck in left and right images
                [nR(n)] = get_neck_edge(obj, xR_filt, right, false);
                [nL(n)] = get_neck_edge(obj, xL_filt, left, false);
                
                if n == 527
                    disp('here')
                end
                
                %min_neck_point = round(1.1*min([nL(n).peak , nR(n).peak]));
                %min_neck_image = obj.head_bw_vid(min_neck_point:obj.head_dim(1),:,n);
                %[~,cent_x] = find(min_neck_image);
                %obj.pivot(n,1) = median(cent_x);
                
                % Get bottom edge points from left & right images
                neck_x = [nL(n).x ; nR(n).x];
                neck_y = [nL(n).y ; nR(n).y];
                obj.neck(n).x = neck_x(~isnan(neck_x));
                obj.neck(n).y = neck_y(~isnan(neck_y));
                
                % Fit a line to the bottom edge points and find the angle in degrees
                c = polyfit(obj.neck(n).x, obj.neck(n).y, 1);
                obj.neck(n).x_test = (1:obj.head_dim(2))';
                obj.neck(n).y_test = polyval(c, obj.neck(n).x_test);
                obj.neck(n).slope = c(1);
                obj.neck(n).yaw = rad2deg(atan( obj.neck(n).slope ));
                
                % Estimate the pivot point of the neck about yaw by finding where the yaw line intersects
                % the image vertical midline
                [~,obj.pivot(n,1)] = min(abs(obj.neck(n).x_test - obj.head_cent(1)));
                %[~,obj.pivot(n,1)] = min(abs(obj.neck(n).x_test - obj.pivot(n,1)));
                obj.pivot(n,2) = obj.neck(n).y_test(obj.pivot(n,1));
                
                %[test,test1] = polyxpoly(obj.neck(n).x_test, obj.neck(n).y_test, ...
                    %[obj.pivot(n,1) obj.pivot(n,1)], [1 obj.head_dim(2)] );
                    
            	% Stabilize the frame in head reference frame
                obj.head_stable_vid(:,:,n) = imrotate(raw, obj.neck(n).yaw, 'bilinear', 'crop');
            end
            obj.neck_left = nL;
            obj.neck_right = nR;
            yaw = cat(1,obj.neck(:).yaw);
            pivot = median(obj.pivot,1);
            obj.yaw = yaw;
            obj.pivot = pivot;
            
            % Estimate when saturation occrurs
            [obj] = detect_saturation(obj);            
            
            fprintf(' Done. \n')
        end
        
        function [neck] = get_neck_edge(obj, edge, I, debug)
            % get_neck_edge:
            %   INPUTS:
            %       edge    : far edge for each row
            %       I       : left or right image
            %       debug   : show plots
            %
            %   OUTPUTS:
            %       neck    : neck structure
            %
            neck.I = I;
            neck.full_edge = edge;
            Isz = size(I);
            [pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                                            'MinPeakProminence', 10, 'SortStr', 'descend');
            idx = (1:length(edge))';
            if isempty(locs)
                %warning('Small peak')
                neck.flag = true;
                [pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                                        'MinPeakProminence', 5);
                if isempty(locs)
                    warning('Can''t find angle')
                    %pause(0.1)
                    neck.flag = nan;
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
            else
                neck.flag = false;
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
                    patch([1 1 Isz(2) Isz(2)], [neck.edge_range(1) neck.edge_range(end) ...
                        neck.edge_range(end) neck.edge_range(1)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                    plot(edge,idx,'m', 'LineWidth', 1)
                    plot(neck.dx, neck.peak, '.g', 'MarkerSize', 15)
                    plot(neck.x, neck.y, '.b', 'LineWidth', 1)
                    plot(neck.x_test, neck.y_test, '-r', 'LineWidth', 1)
                set(ax, 'Box', 'on')
            end
        end
        
        function [obj] = detect_saturation(obj)
            % detect_saturation: find roll angle of head in each frame
            %  
            %
            
            left = cat(1, obj.neck_left(:).dx);
            right = cat(1, obj.neck_right(:).dx);
            [obj.saturation(:,1)] = isoutlier(left, 'quartiles', 1);
            [obj.saturation(:,2)] = isoutlier(right, 'quartiles', 1);
        end
        
        function [roll, obj] = get_roll(obj, eye_ratio, roll_cal)
            % get_roll: find roll angle of head in each frame
            %  eye_ratio: amount around centroid to get eye intensity
            %  roll_cal: roll calibration factor
            %
            
            obj.roll_cal = roll_cal;
            
            % Isolate head in each frame
            [obj] = iso_head(obj, 2);
            
            % Get head radiius
           	bb = cat(1,obj.head_iso_props(:).BoundingBox);
            obj.head_radius = median(bb(:,4));
            obj.head_tip = obj.pivot + obj.head_radius * [sind(obj.yaw) , -cosd(obj.yaw)];
            
            % Find window to measure eye intensity
           	obj.eye.cent = cat(1,obj.head_iso_props(:).Centroid);
            obj.eye.cent_med = round(median(obj.eye.cent,1));
            obj.eye.span = round(eye_ratio * obj.head_iso_props(1).BoundingBox(4) / 2);
            obj.eye.window = (obj.eye.cent_med(2) - obj.eye.span : obj.eye.cent_med(2) + obj.eye.span)';
            
            % Find roll angles
            fprintf('Calculating roll angles ...')
            obj.eye.intensity = cell(obj.head_dim(3),1);
            obj.eye.intensity_mean = nan(obj.head_dim(3),obj.head_dim(2));
            for n = 1:obj.head_dim(3)
                %fprintf('%i \n', n)
                frame = obj.head_iso_vid(:,:,n); % frame for roll analysis
                
                % Get eye intensityies in window
             	obj.eye.intensity{n} = double(frame(obj.eye.window,:)); % intensity in eye range
                obj.eye.intensity_mean(n,:) = mean(obj.eye.intensity{n},1); % mean intensity
                
                % Find peaks
                [pks,locs,w,p] = findpeaks(obj.eye.intensity_mean(n,:), 'MinPeakHeight', 100, 'MinPeakDistance', 5, ...
                                'MinPeakWidth', 5, 'MinPeakProminence', 5);
                            
                % Get left and right eye peaks (1st & last)
            	obj.eye.peak_loc(n,:) = locs([1,end]); % pixel in horizontal plane
                obj.eye.peak_int(n,:) = pks([1,end]); % intensity value (close too 255)
                obj.eye.search_int(n,1) = min(0.95*obj.eye.peak_int(n,:)); % intensity level to find eye widths
                
                % Find where left and right curves intersect the search level
                [xi,~] = polyxpoly([1 obj.head_dim(2)], [obj.eye.search_int(n,1) obj.eye.search_int(n,1)], ...
                    1:obj.head_dim(2), obj.eye.intensity_mean(n,:));
                obj.eye.left(n,:) = [xi(1) , xi(2)];
                obj.eye.right(n,:) = [xi(end-1) , xi(end)];
            end
            % Find eye widths and roll index
          	obj.eye.left_width  = diff(obj.eye.left, 1, 2);
            obj.eye.right_width = diff(obj.eye.right, 1, 2);
         	obj.roll_idx = (obj.eye.left_width - obj.eye.right_width) ./ ...
                (obj.eye.left_width + obj.eye.right_width);
            roll = obj.roll_idx .* obj.roll_cal;
            obj.roll = roll;
            fprintf(' Done. \n')
        end   
        
        function [obj] = iso_head(obj, med_win)
            % iso_head: isolate head in yaw stabilized video
            %  med_win: median filter window size
            %
            
            fprintf('Isolating head ...')
            
            % Filter video and get find where the head starts in every frame
            vid_filt = obj.head_stable_vid;
            obj.head_start = nan(obj.head_dim(3),1);
            for n = 1:obj.head_dim(3)
                frame = obj.head_stable_vid(:,:,n);
                vid_filt(:,:,n) = medfilt2(frame, med_win*[1 1]);
                frame_bw = imbinarize(vid_filt(:,:,n));
                obj.head_start(n) = find( max(frame_bw,[],2) , 1, 'first');
            end
            obj.head_start_min = min(obj.head_start);
            
            % Isolate head in each frame
            y_range = obj.head_start_min:obj.pivot(2);
            obj.head_iso_vid = uint8( zeros( length(y_range), obj.head_dim(2)) );
          	for n = 1:obj.head_dim(3)
                obj.head_iso_vid(:,:,n) = vid_filt(y_range,:,n);
                
                bw_iso = imbinarize(obj.head_iso_vid(:,:,n));
                image_props = regionprops(bw_iso,'Centroid','BoundingBox','MajorAxisLength');
                
                if length(image_props) > 1
                    mjaxes = [image_props.MajorAxisLength];
                    [~,idx] = sort(mjaxes, 'descend');
                    idx = idx(1);
                else
                    idx = 1;
                end  
                
                vid_props_iso(n) = image_props(idx);
            end
            obj.head_iso_props = vid_props_iso;
            fprintf(' Done. \n')
        end
        
        function [obj] = play_tracking(obj, playback, tt, targetdir, fname)
            % play_tracking: show tracking anotations
            %   playback: playback rate
            %   tt: time vector
            %   targetdir: target directory for output video
            %   fname: file name
            %

           	if ~playback
                return
            end
            
            if nargin < 4
                export = false;
                if nargin < 3
                    tt = (1:obj.head_dim(3))';
                end
            else
                export = true;
                [~,outname,ext] = fileparts(fname);
                outpath = fullfile(targetdir, [outname '_montage_head.' ext]);
                if ~isempty(tt)
                    Fs = round( 1 / mean(diff(tt)));
                    vidFs = round( Fs / playback);
                else
                    tt = (1:obj.head_dim(3))';
                    vidFs = 50;
                    warning('No time given, using 50 fps')
                end
                VID = VideoWriter(outpath,'MPEG-4');
                VID.FrameRate = vidFs;
                open(VID)
            end
            
        	fig = figure ; clf
            set(fig, 'Color', [0 0 0.1], 'Units', 'inches', 'Position', [4 1 10 9])
                ax(1) = subplot(4,3,2); cla ; hold on ; axis image
                ax(2) = subplot(4,3,5); cla ; hold on ; axis image
                ax(3) = subplot(4,3,1); cla ; hold on ; axis image
                ax(4) = subplot(4,3,3); cla ; hold on ; axis image
                ax(5) = subplot(4,3,4); cla ; hold on ; axis image
                ax(6) = subplot(4,3,6); cla ; hold on ; axis image
                ax(7) = subplot(4,3,7); cla ; hold on ; axis image
                ax(8) = subplot(4,3,8); cla ; hold on ; axis image
                ax(9) = subplot(4,3,9); cla ; hold on ; axis tight ; ylabel('Intensity') ; ylim([100 270])
                ax(10) = subplot(4,3,10:12); cla ; hold on ; ylabel('Angle (°)') ; xlabel('Frame')
                    h.yaw = animatedline(ax(10), 'Color', 'c', 'LineWidth', 1);
                    h.roll = animatedline(ax(10), 'Color', 'r', 'LineWidth', 1);
                    xlim([0 tt(end)])
                    
                set(ax(1:10), 'XColor', 'w', 'YColor', 'w', 'Box', 'on', 'FontSize', 9)
                set(ax(1:8), 'Color', 'none')
                set(ax(9), 'Color', 'k')
                set(ax(10), 'Color', 'k', 'LineWidth', 1, 'Box', 'off')
                set(ax(1:9), 'XTick', [1 20:20:obj.head_dim(2)], 'XLim', [0 obj.head_dim(2)])
                set(ax(1:8), 'YTick', [1 20:20:obj.head_dim(1)])
                set(ax(1:9), 'Visible', 'on')
                
            for n = 1:obj.head_dim(3)
                ax(10) = subplot(4,3,10:12); hold on
                    addpoints(h.yaw, tt(n), obj.yaw(n))
                    addpoints(h.roll, tt(n), obj.roll(n))
                    drawnow
                    pause(0.0001)
                    %set(ax(1:9), 'Visible', 'off')
                    
                if ~mod(n+1,playback)
                    ax(1) = subplot(4,3,2); cla ; hold on
                    str = sprintf('Yaw = %2.3f°', obj.yaw(n));
                    title(str, 'Color', 'c', 'FontWeight', 'bold', 'FontSize', 14)                    
                        imshow(obj.head_vid(:,:,n))
                        plot([obj.pivot(1) obj.pivot(1)], [1 obj.head_dim(1)], '--m', 'LineWidth', 0.5)
                        plot([1 obj.head_dim(2)], [obj.pivot(2) obj.pivot(2)], '--m', 'LineWidth', 0.5)
                        plot(obj.neck(n).x, obj.neck(n).y, '.b', 'MarkerSize', 5)
                        plot(obj.neck(n).x_test, obj.neck(n).y_test, 'c', 'LineWidth', 2)
                        plot([obj.pivot(1) obj.head_tip(n,1)], [obj.pivot(2) obj.head_tip(n,2)], ...
                            'c', 'LineWidth', 2)
                        plot(obj.pivot(1), obj.pivot(2), '.c', 'MarkerSize', 15)
                        
                    ax(2) = subplot(4,3,5); cla ; hold on
                        imshow(obj.head_out_vid(:,:,n))
                        plot([obj.pivot(1) obj.pivot(1)], [1 obj.head_dim(1)], '--m', 'LineWidth', 0.5)
                        plot([1 obj.head_dim(2)], [obj.pivot(2) obj.pivot(2)], '--m', 'LineWidth', 0.5)
                        plot(obj.neck(n).x, obj.neck(n).y, '.b', 'MarkerSize', 5)
                        plot(obj.neck(n).x_test, obj.neck(n).y_test, 'c', 'LineWidth', 1)
                        plot([obj.pivot(1) obj.head_tip(n,1)], [obj.pivot(2) obj.head_tip(n,2)], ...
                            'c', 'LineWidth', 1)
                        plot(obj.pivot(1), obj.pivot(2), '.c', 'MarkerSize', 15)
                        plot(obj.pivot(1), obj.neck_left(n).peak, '.g', 'MarkerSize', 15)
                      	plot(obj.pivot(1), obj.neck_right(n).peak, '.y', 'MarkerSize', 15)

                    ax(3) = subplot(4,3,1); cla ; hold on
                        imshow(0*obj.head_out_vid(:,:,n))
                        if obj.saturation(n,1)
                            title('Contact!', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'r')
                            patch([1 1 obj.head_dim(2) obj.head_dim(2)], ...
                                [1 obj.head_dim(1) obj.head_dim(1) 1], 'r', ...
                                'FaceAlpha', 0.2, 'EdgeColor', 'none')
                        else
                            title('')
                        end
                        plot(flipud(obj.neck_left(n).full_edge), obj.head_dim(1):-1:1, ...
                            'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
                        plot(obj.neck_left(n).dx, obj.neck_left(n).peak, '.g', 'MarkerSize', 15)

                    ax(4) = subplot(4,3,3); cla ; hold on
                        imshow(0*obj.head_out_vid(:,:,n))
                        if obj.saturation(n,2)
                            title('Contact!', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'r')         
                            patch([1 1 obj.head_dim(2) obj.head_dim(2)], ...
                                [1 obj.head_dim(1) obj.head_dim(1) 1], 'r', ...
                                'FaceAlpha', 0.2, 'EdgeColor', 'none')
                        else
                            title('')
                        end
                        plot(flipud(obj.neck_right(n).full_edge), obj.head_dim(1):-1:1, ...
                            'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
                        plot(obj.neck_right(n).dx, obj.neck_right(n).peak, '.y', 'MarkerSize', 15)
                        
                    ax(5) = subplot(4,3,4); cla ; hold on
                        imshow(obj.neck_left(n).I)
                        patch([1 1 obj.head_dim(2) obj.head_dim(2)], [obj.neck_left(n).edge_range(1) obj.neck_left(n).edge_range(end) ...
                            obj.neck_left(n).edge_range(end) obj.neck_left(n).edge_range(1)], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                        plot(obj.neck_left(n).x, obj.neck_left(n).y, '.b', 'MarkerSize', 5)
                        %plot(obj.neck_left(n).x_test, obj.neck_left(n).y_test, 'y', 'LineWidth', 1)

                    ax(6) = subplot(4,3,6); cla ; hold on
                        imshow(obj.neck_right(n).I)
                        patch([1 1 obj.head_dim(2) obj.head_dim(2)], [obj.neck_right(n).edge_range(1) obj.neck_right(n).edge_range(end) ...
                            obj.neck_right(n).edge_range(end) obj.neck_right(n).edge_range(1)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                        plot(obj.neck_right(n).x, obj.neck_right(n).y, '.b', 'MarkerSize', 5)
                        %plot(obj.neck_right(n).x_test, obj.neck_right(n).y_test, 'y', 'LineWidth', 1)
                        
                    ax(7) = subplot(4,3,7); cla ; hold on
                        imshow(obj.head_stable_vid(:,:,n))
                        plot([1 obj.head_dim(2)], [obj.pivot(2) obj.pivot(2)], '--c', 'LineWidth', 0.5)
                        plot([1 obj.head_dim(2)], [obj.head_start_min obj.head_start_min], '--c', 'LineWidth', 0.5)
                        plot(obj.pivot(1), obj.pivot(2), '.c', 'MarkerSize', 15)

                    ax(8) = subplot(4,3,8); cla ; hold on
                    str = sprintf('Roll = %2.3f°', obj.roll(n));
                    title(str, 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 14)  
                        imshow(obj.head_iso_vid(:,:,n))
                        x = [1 1 obj.head_dim(2) obj.head_dim(2)];
                        y = [obj.eye.window(1) obj.eye.window(end) obj.eye.window(end) obj.eye.window(1)];
                        patch(x, y, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
                        plot(obj.eye.left(n,:), [obj.eye.cent_med(2) obj.eye.cent_med(2)],  '-r', 'LineWidth', 3)
                        plot(obj.eye.right(n,:), [obj.eye.cent_med(2) obj.eye.cent_med(2)],  '-r', 'LineWidth', 3)
                        
                    ax(9) = subplot(4,3,9); cla ; hold on
                        plot(obj.eye.intensity{n}', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
                        plot(obj.eye.intensity_mean(n,:), 'w', 'LineWidth', 2)
                        plot(obj.eye.peak_loc(n,1), obj.eye.peak_int(n,1), '.g', 'MarkerSize', 20)
                        plot(obj.eye.peak_loc(n,2), obj.eye.peak_int(n,2), '.g', 'MarkerSize', 20)
                        plot([1 obj.head_dim(2)], [obj.eye.search_int(n) obj.eye.search_int(n)], '-b', 'LineWidth', 1)
                        plot(obj.eye.left(n,:), [obj.eye.search_int(n) obj.eye.search_int(n)],  '-r', 'LineWidth', 3)
                        plot(obj.eye.right(n,:), [obj.eye.search_int(n) obj.eye.search_int(n)],  '-r', 'LineWidth', 3) 
                
                    if export
                        fig_frame = getframe(fig);
                        writeVideo(VID,fig_frame);
                    end
                end
                pause(0.0001)
            end
        end
        
    end
    
end











