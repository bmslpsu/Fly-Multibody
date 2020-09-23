classdef headtracker_area_v2
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
        yaw_filt            % filtered yaw for head stabilization (to compute roll)
        saturation          % saturation [L, R] (boolean)
        
        eye                 % eye properties
        roll_idx            % head roll ratio
        roll_cal            % head roll calibration factor
        roll                % head roll
    end
    
    methods
        function obj = headtracker_area_v2(vid)
            % headtracker_area: Construct an instance of this class
            %   
            
            obj.vid = squeeze(vid); % input video
            obj.dim = size(obj.vid); % video size
            
            % Get video properties
            [obj.vid_props, obj.bw_vid, obj.bound_box, obj.body_yaw] = get_vid_props(obj, vid, [], false);
            
            % Get head reigon
            [obj] = get_head_vid(obj, [0 0.1], 0.3);
            
            % Get yaw
            [~, ~, obj] = get_yaw(obj, 0.5);
            
            % Isolate head in each frame
            [obj] = iso_head(obj, 2);
            
            % Get roll
            [~, obj] = get_roll(obj, 0.35, 36.33, false);
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
            se_close = strel('disk',1);
            for n = 1:obj.dim(3)
                A = cell(1);
                A{1} = 2*obj.head_vid(:,:,n);
                %A{end+1} = medfilt2(A{end}, [5 5]);
                A{end+1} = imbinarize(A{end});
                A{end}(end,:) = true;
                A{end+1} = imclose(A{end}, se_close);
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
                A{end+1} = bwmorph(A{end}, 'bridge',3);
                
                %imshow(A{end})
                
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
                %raw = obj.head_vid(:,:,n); % raw for stabilization
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
                
                if n == 120
                    disp('here')
                end
                
                % Get bottom edge of neck in left and right images
                if n > 1
                    [nR(n)] = get_neck_edge(obj, xR_filt, right, nR, false);
                    [nL(n)] = get_neck_edge(obj, xL_filt, left, nL, false);
                else
                    [nR(n)] = get_neck_edge(obj, xR_filt, right, [], false);
                    [nL(n)] = get_neck_edge(obj, xL_filt, left, [], false);
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
            end
            obj.neck_left = nL;
            obj.neck_right = nR;
            
            yaw = cat(1,obj.neck(:).yaw);
            %yaw = hampel(1:obj.head_dim(3), yaw);
            obj.yaw = yaw;
            
            pivot = median(obj.pivot,1);
            obj.pivot = pivot;
            
            % Estimate when saturation occrurs
            [obj] = detect_saturation(obj);            
            
            fprintf(' Done. \n')
        end
        
        function [neck] = get_neck_edge(obj, edge, I, side, debug)
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
            cut_edge = edge;
            cutI = round(0.3 * obj.head_dim(1));
            cut_edge(1:cutI) = nan;
            Isz = size(I);
            [pks,locs,w,p] = findpeaks(-cut_edge, 'MinPeakDistance', 10, 'MinPeakWidth', 3, ...
                                            'MinPeakProminence', 10, 'SortStr', 'descend');
            idx = (1:length(edge))';
            if isempty(locs)
                %warning('Small peak')
                neck.flag = true;

                if ~isempty(side)
                    good_peak = ~[side.flag]';
                    pks = -[side.dx];
                    pks = median(pks(good_peak));
                    locs = [side.peak];
                    locs = median(locs(good_peak));
                    w = [side.width];
                    w = median(w(good_peak));
                else
                    [pks,locs,w,p] = findpeaks(-edge, 'MinPeakDistance', 10, ...
                            'MinPeakWidth', 3, 'MinPeakProminence', 5);
                end
                
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
            
            neck.peak = round(locs(1));
            neck.dx = -pks(1);
            neck.width = w(1);
            
            if ~neck.flag
                offset = 0.2*neck.width;
                Lr = 1.5;
            else
                offset = -0.4*neck.width;
                Lr = 2;
            end

            neck.edge_range = ((neck.peak  - ceil(neck.width/Lr) ):floor(neck.peak - offset))';
            neck.edge = (edge(neck.edge_range));
            neck.image = I(neck.edge_range,:);
            neck.image_all = false(size(I));
            neck.image_all(neck.edge_range,:) = neck.image;

            [y, x] = find(neck.image_all);
            neck.x = x;
            neck.y = y;

            if debug
                fig = figure (109); cla
                set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 4 6])
                ax(1) = subplot(1,1,1); hold on; axis image
                    imshow(255*I, 'InitialMagnification', 400) ; hold on
                    patch([1 1 Isz(2) Isz(2)], [neck.edge_range(1) neck.edge_range(end) ...
                        neck.edge_range(end) neck.edge_range(1)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                    plot(edge,idx,'m', 'LineWidth', 1)
                    plot(neck.dx, neck.peak, '.g', 'MarkerSize', 15)
                    plot(neck.x, neck.y, '.b', 'LineWidth', 1)
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
        
        function [roll, obj] = get_roll(obj, eye_ratio, roll_cal, debug)
            % get_roll: find roll angle of head in each frame
            %  eye_ratio: amount around centroid to get eye intensity
            %  roll_cal: roll calibration factor
            %
            res = 5;
            obj.roll_cal = roll_cal;
                        
            % Get head radiius
           	bb = cat(1,obj.head_iso_props(:).BoundingBox);
            obj.head_radius = median(bb(:,4));
            obj.head_tip = obj.pivot + obj.head_radius * [sind(obj.yaw) , -cosd(obj.yaw)];
            
            % Find window to measure eye intensity
            obj.eye = [];
            obj.eye.xspan = 1:obj.head_dim(2);
            obj.eye.xinterp = 1:1/res:obj.head_dim(2);
           	obj.eye.cent = cat(1,obj.head_iso_props(:).Centroid);
            obj.eye.cent_med = round(median(obj.eye.cent,1));
            obj.eye.span = round(eye_ratio * obj.head_iso_props(1).BoundingBox(4) / 2);
            obj.eye.window = (obj.eye.cent_med(2) - obj.eye.span : obj.eye.cent_med(2) + obj.eye.span)';
            
            % Find roll angles
            fprintf('Calculating roll angles ...')
            obj.eye.intensity = cell(obj.head_dim(3),1);
            %obj.eye.intensity_mean = nan(obj.head_dim(3),obj.head_dim(2));
            [b,a] = butter(3, 0.3, 'low');
            for n = 1:obj.head_dim(3)
                %fprintf('%i \n', n)
                frame = obj.head_iso_vid(:,:,n); % frame for roll analysis
                %frame = interp_image(obj, frame, 10);
                
                % Get eye intensities in window
             	obj.eye.intensity{n} = double(frame(obj.eye.window,:)); % intensity in eye range
                int_mean = filtfilt(b, a, mean(obj.eye.intensity{n},1));
                int_mean_interp = interp1(obj.eye.xspan, int_mean, obj.eye.xinterp, 'pchip');
                obj.eye.intensity_mean(n,:) = int_mean_interp; % mean intensity
                dx = diff(obj.eye.intensity_mean(n,:)) / mean(diff(obj.eye.xinterp));
                dx = [dx(1) , dx];
                dx = filtfilt(b, a, dx);
                obj.eye.intensity_mean_diff(n,:) = dx; % derivative of mean intensity

                % Find outer spikes to get eye starts (outer)
                [~,locs,w] = findpeaks(obj.eye.intensity_mean_diff(n,:), 'MinPeakHeight', 2, ...
                                'MinPeakWidth', res*2, 'MinPeakProminence', 0.5);
              	left_outer_I = (locs(1) - 0*w(1)/1.5);
                left_outer_x = obj.eye.xinterp(left_outer_I);
                
                [~,locs,w] = findpeaks(-obj.eye.intensity_mean_diff(n,:), 'MinPeakHeight', 2, ...
                                'MinPeakWidth', res*2, 'MinPeakProminence', 0.5);
               	right_outer_I = (locs(end) + 0*w(end)/1.5);
                right_outer_x = obj.eye.xinterp(right_outer_I);
                
                % Find peak intensities
                [pks,locs,w,p] = findpeaks(obj.eye.intensity_mean(n,:), 'MinPeakHeight', 100, ...
                                'MinPeakDistance', res*5, 'MinPeakWidth', res*2, 'MinPeakProminence', 5);

                % Get left and right eye peaks (1st & last)
            	obj.eye.peak_loc(n,:) = locs([1,end]); % pixel in horizontal plane
                obj.eye.peak_x(n,:) = obj.eye.xinterp(locs([1,end])); % interp pixel in horizontal plane
                obj.eye.peak_int(n,:) = pks([1,end]); % intensity value (close too 255)
                
                % Get change in intensity between peaks
                flex = round(0.05 * obj.head_dim(2));
                obj.eye.eye_range{n} = obj.eye.peak_loc(n,1)-flex : obj.eye.peak_loc(n,2)+flex; % between peaks
                mid = nan(size(dx));
                mid(obj.eye.eye_range{n}) = obj.eye.intensity_mean_diff(n,obj.eye.eye_range{n});
                obj.eye.eye_mid{n} = mid; % intensity between peaks
                
              	[~,locs,w] = findpeaks(-obj.eye.eye_mid{n}, 'MinPeakHeight', 1, ...
                                'MinPeakWidth', res*2, 'MinPeakProminence', 3);
              	if ~isempty(locs)
                    left_inner = (locs(1) + 0*w(1)/1.5);
                    obj.eye.left_inner_int(n) = obj.eye.intensity_mean(n,round(left_inner));
                    obj.eye.left(n,:) = [left_outer_x, obj.eye.xinterp(left_inner)];
                else
                    obj.eye.left_inner_int(n) = 0;
                    obj.eye.left(n,:) = [0, 0];
                end
                
             	[~,locs,w] = findpeaks(obj.eye.eye_mid{n},'MinPeakHeight', 1, ...
                                'MinPeakWidth', res*2, 'MinPeakProminence', 3);
               	if ~isempty(locs)
                    right_inner = (locs(end) - 0*w(end)/1.5);
                    obj.eye.right_inner_int(n) = obj.eye.intensity_mean(n,round(right_inner));
                    obj.eye.right(n,:) = [obj.eye.xinterp(right_inner) , right_outer_x];
                else
                    obj.eye.right_inner_int(n) = 0;
                    obj.eye.right(n,:) = [0 , 0];
                end
                
                if debug
                    figure (111) ; clf ; hold on ; ylim([-50 255])
                        plot(obj.eye.xinterp, obj.eye.intensity_mean(n,:), 'k', 'LineWidth', 1)
                        plot(obj.eye.xinterp([1 end]), [0 0], '--k', 'LineWidth', 0.5)
                        plot(obj.eye.peak_x(n,:), obj.eye.peak_int(n,:), '.r', 'MarkerSize', 15)
                        plot(obj.eye.xinterp, obj.eye.intensity_mean_diff(n,:), 'g', 'LineWidth', 1)
                        plot(obj.eye.xinterp, obj.eye.eye_mid{n}, 'r', 'LineWidth', 1.5)
                        plot(obj.eye.left(n,:), [obj.eye.left_inner_int(n) obj.eye.left_inner_int(n)], ...
                            '.-m', 'MarkerSize', 15)
                        plot(obj.eye.right(n,:), [obj.eye.right_inner_int(n) obj.eye.right_inner_int(n)], ...
                            '.-c', 'MarkerSize', 15)
                        plot([obj.eye.left(n,1) obj.eye.left(n,1)], [0 obj.eye.left_inner_int(n)], 'm')
                        plot([obj.eye.left(n,2) obj.eye.left(n,2)], [0 obj.eye.left_inner_int(n)], 'm')
                        plot([obj.eye.right(n,2) obj.eye.right(n,2)], [0 obj.eye.right_inner_int(n)], 'c')
                        plot([obj.eye.right(n,1) obj.eye.right(n,1)], [0 obj.eye.right_inner_int(n)], 'c')
                end
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
            [b,a] = butter(3, 0.3, 'low');
            obj.yaw_filt = obj.yaw;
            obj.yaw_filt = hampel(1:obj.dim(3), obj.yaw_filt, 10, 3, 'Adaptive', 0.1);
            obj.yaw_filt = filtfilt(b,a,obj.yaw_filt);
            
           	% Stabilize the frame in head reference
            for n = 1:obj.head_dim(3)
                obj.head_stable_vid(:,:,n) = imrotate(obj.head_vid(:,:,n), ...
                    obj.neck(n).yaw, 'bilinear', 'crop');
            end
            
            % Filter video and find where the head starts in every frame
            vid_filt = obj.head_stable_vid;
            obj.head_start = nan(obj.head_dim(3),1);
            for n = 1:obj.head_dim(3)
                frame = obj.head_stable_vid(:,:,n);
                vid_filt(:,:,n) = medfilt2(frame, med_win*[1 1]);
                frame_bw = imbinarize(vid_filt(:,:,n));
                obj.head_start(n) = find( max(frame_bw,[],2) , 1, 'first');
            end
            obj.head_start_min = median(obj.head_start);
            
            % Isolate head in each frame
            y_range = obj.head_start_min:obj.pivot(2);
            obj.head_iso_vid = uint8( zeros( length(y_range), obj.head_dim(2)) );
          	for n = 1:obj.head_dim(3)
                obj.head_iso_vid(:,:,n) = obj.head_stable_vid(y_range,:,n);
                
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
        
       	function [I_interp] = interp_image(obj, I, res)
            % interp_image: interpolate a 2D image
            %   I: input image
            %   res: amount to scale down pixel (1/res)
            %
            
            step = 1 / res;
            dimI = size(I);
            class_of_I = class(I);
            [x, y] = meshgrid(1:dimI(2), 1:dimI(1));
            [xi, yi] = meshgrid(1:step:dimI(2), 1:step:dimI(1));
            I_interp = cast( interp2(x, y, double(I), xi, yi,'linear'), class_of_I);
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
                outpath = fullfile(targetdir, [outname '_montage_head' ext]);
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
            set(fig, 'Color', [0 0 0], 'Units', 'inches', 'Position', [4 1 10 9])
                ax(1) = subplot(4,3,2); cla ; hold on ; axis image
                ax(2) = subplot(4,3,5); cla ; hold on ; axis image
                ax(3) = subplot(4,3,1); cla ; hold on ; axis image
                ax(4) = subplot(4,3,3); cla ; hold on ; axis image
                ax(5) = subplot(4,3,4); cla ; hold on ; axis image
                ax(6) = subplot(4,3,6); cla ; hold on ; axis image
                ax(7) = subplot(4,3,7); cla ; hold on ; axis image
                ax(8) = subplot(4,3,8); cla ; hold on ; axis image
                ax(9) = subplot(4,3,9); cla ; hold on ; axis tight ; ylabel('Intensity') ; xlabel('Frame Width') ; ylim([-50 255])
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
%                             title('Contact!', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'r')
%                             patch([1 1 obj.head_dim(2) obj.head_dim(2)], ...
%                                 [1 obj.head_dim(1) obj.head_dim(1) 1], 'r', ...
%                                 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                        else
                            title('')
                        end
                        plot(flipud(obj.neck_left(n).full_edge), obj.head_dim(1):-1:1, ...
                            'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
                        plot(obj.neck_left(n).dx, obj.neck_left(n).peak, '.g', 'MarkerSize', 15)

                    ax(4) = subplot(4,3,3); cla ; hold on
                        imshow(0*obj.head_out_vid(:,:,n))
                        if obj.saturation(n,2)
%                             title('Contact!', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'r')         
%                             patch([1 1 obj.head_dim(2) obj.head_dim(2)], ...
%                                 [1 obj.head_dim(1) obj.head_dim(1) 1], 'r', ...
%                                 'FaceAlpha', 0.2, 'EdgeColor', 'none')
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
                        plot(obj.eye.xspan, obj.eye.intensity{n}', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
                        plot(obj.eye.xinterp, obj.eye.intensity_mean(n,:), 'w', 'LineWidth', 2)
                        plot(obj.eye.xinterp([1 end]), [0 0], '--y', 'LineWidth', 0.5)
                        plot(obj.eye.peak_x(n,:), obj.eye.peak_int(n,:), '.g', 'MarkerSize', 15)
                        plot(obj.eye.xinterp, obj.eye.intensity_mean_diff(n,:), 'b', 'LineWidth', 1)
                        plot(obj.eye.xinterp, obj.eye.eye_mid{n}, 'g', 'LineWidth', 1.5)
                        plot(obj.eye.left(n,:), [obj.eye.left_inner_int(n) obj.eye.left_inner_int(n)], ...
                            '.-r', 'MarkerSize', 15)
                        plot(obj.eye.right(n,:), [obj.eye.right_inner_int(n) obj.eye.right_inner_int(n)], ...
                            '.-r', 'MarkerSize', 15)                  
                     	 plot([obj.eye.left(n,1) obj.eye.left(n,1)], [0 obj.eye.left_inner_int(n)], 'r')
                        plot([obj.eye.left(n,2) obj.eye.left(n,2)], [0 obj.eye.left_inner_int(n)], 'r')
                        plot([obj.eye.right(n,2) obj.eye.right(n,2)], [0 obj.eye.right_inner_int(n)], 'r')
                        plot([obj.eye.right(n,1) obj.eye.right(n,1)], [0 obj.eye.right_inner_int(n)], 'r')
                        %plot([1 obj.head_dim(2)], [obj.eye.search_int(n) obj.eye.search_int(n)], '-b', 'LineWidth', 1)
                        %plot(obj.eye.left(n,:), [obj.eye.search_int(n) obj.eye.search_int(n)],  '-r', 'LineWidth', 3)
                        %plot(obj.eye.right(n,:), [obj.eye.search_int(n) obj.eye.search_int(n)],  '-r', 'LineWidth', 3) 
                
                    if export
                        fig_frame = getframe(fig);
                        writeVideo(VID,fig_frame);
                    end
                end
                pause(0.0001)
            end
        end
        
    end
    
    methods (Static)
        
    end
    
end











