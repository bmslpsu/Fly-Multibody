classdef saccade_v1
    % saccade_v1: class to detect and analyze saccades in signal
    %   
 	
    properties (SetAccess = private, Hidden = false)
        % Input signal attributes
        time                    % raw time
        position                % raw position
        velocity                % raw velocity
        acceleration            % raw acceleration
        abs_velocity            % absolute value of velocity
        position_filt_detect    % low- & high-pass filtered position for saccade peak detection
        velocity_filt_detect    % low- & high-pass filtered velocity for saccade peak detection
        position_filt_ss        % low- & high-pass filtered position for saccade start/stop detection
        velocity_filt_ss        % low- & high-pass filtered velocity for saccade start/stop detection
     	n                       % # of data points
        Ts                      % sampling time
        Fs                      % sampling frequency
     	index                   % indicies
        
        % Saccade detection
     	ismanual                % boolean: true if saccades were detetected manually and passed as "peaks"
        direction               % direction of saccades to detect (- = -1, all = 0, + = 1)
      	amp_cut                 % cut off for saccade amplitudes (< amp_cut are removed)
        dur_cut                 % cut off for saccade durations (> dur_cut are removed)
        sacd_length             % make every saccade the same length if not nan
    	Fc_detect               % 2x1 low- & high-pass filter frequencies for detecing saccade peaks
        Fc_ss                   % 2x1 low- & high-pass filter frequencies for detecing saccade starts/ends
     	nstd_thresh             % # of std's used for detection threshold (nan if manual)
        threshold               % saccade detection threshold (nan if manual)   
       	true_thresh             % cut-off for peak velocity in unfiltered raw signal
        pks                     % saccade peak locations [indicies] (optional)
        min_pkdist              % minimum peak distance [s]
        min_pkwidth             % minimum peak width [s]
        min_pkprom              % minimum peak promiance
        min_pkthresh            % minimum height difference between peak and its neighbors
        boundThresh             % ratio of saccade peak velocties to determine start and stop
        
        % Saccade atrributes
        count                   % # of saccades deteteced
        rate                    % saccades / time unit
        sacd_ratio              % saccade time / total time
        SACD                    % table of saccade statistics
        
        % Isolated saccades and inter-saccade intervals
        saccades                % saccade only index/time/position/velocity tables in cells
        intervals               % inter-saccade intervals index/time/position/velocity tables in cells
        shift                   % data with interpolation between saccades
       	saccades_all            % all saccade data
        removed_all             % data with saccades removed & nan's in their place
       	peaks
        starts
        ends
        
        % Normalized saccades and inter-saccade intervals
        normpeak_saccade        % normalized saccades , aligned to peak time
        norm_interval           % normalized intervals, aligned to start time
        normstart_interval      % normalized intervals, aligned to start time and initial position
        normend_interval        % normalized intervals, aligned to end time and end position

        % Isolated stimulus and error properties during saccades and inter-saccade intervals
        stimlus_position        % stimulus position
        stimlus_velocity        % stimulus veloicty
        stimulus                % stimulus inter-saccade intervals
        normstart_stimulus      % normalized stimulus intervals, aligned to start time
        normend_stimulus        % normalized stimulus intervals, aligned to end time
        error                   % error (stimulus - position) for inter-saccade intervals
        int_error               % integrated error (stimulus - position) for inter-saccade intervals
    end
    
  	properties (SetAccess = public, Hidden = false)
        cmap                    % colormap
        extra                   % anything else to store
    end
    
 	properties (SetAccess = private, Hidden = false)
    	medians
        stds
        amplitudes
        durations
        directions
        tablenames
        
        align_saccade_peak
        align_interval_start
        align_interval_end
    end
     
	properties (Transient)
        
    end
    
    methods
        function obj = saccade_v1(position, time, threshold, true_thresh, Fc_detect, Fc_ss, amp_cut, ...
                                dur_cut, direction, pks, sacd_length, min_pkdist, min_pkwidth, min_pkprom, ...
                                min_pkthresh, boundThresh, showplot)
            % saccade_v1: Construct an instance of this class
            %
            %   INPUTS:
            %       position        : signal vector
            %       time            : time vector
            %       threshold       : absolute velocity threshold for peak detetcion in filtered signal
            %                         (if negative, interpret as STD's from median). If 2x1, use 1st 
            %                         value as STDs from median & 2nd value as minimum threshold.
            %       true_thresh     : cut-off for peak velocity in raw signal
            %       Fc_detect       : 2x1 low- & high-pass filter frequencies for detecing saccade peaks
            %       Fc_ss           : 2x1 low- & high-pass filter frequencies for detecing saccade 
            %                      	  starts and ends
            %       amp_cut         : minimum position amplitude for a saccade
            %       dur_cut         : remove saccades over this length of time
            %       pks             : saccade peak locations [indicies] (optional)
            %       direction       : get saccades in this direction only (1=+, -1=-, 0=both)
            %       sacd_length     : if specified, return each saccade with a set start and stop 
            %                         time offset from peak time
            %       min_pkdist      : minimum peak distance [s]
            %       min_pkwidth     : minimum peak width [s]
            %       min_pkprom      : minimum peak promiance
            %       min_pkthresh    : minimum height difference between peak and its neighbors
            %       boundThresh     : ratio of peaks to determine start and stop
            %
            %       showplot        : show some plots if true
            
            %
            
            if nargin == 0
                obj = setDeafults(obj);
                return
            end
            
            struct_flag = false;
            
            if nargin < 12
                showplot = false; % default: show plots
                if nargin < 11
                    sacd_length = nan; % defualt: find start/stop of saccades automatically
                    if nargin < 10
                        pks = []; % default: no peaks given
                        if nargin < 9
                            direction = 0; % default: both directions
                            if nargin < 8
                                dur_cut = inf;
                                if nargin < 7
                                    amp_cut = 0; % defualt: don't use cutoff
                                    if nargin < 6
                                        Fc_ss = nan(1,2); % defualt: no filter
                                        if nargin < 5
                                             if isstruct(threshold) % properties set by structure
                                               struct_flag = true;
                                               struct_inputs = threshold;
                                               showplot = true_thresh;
                                             else
                                                Fc_detect = nan(1,2); % defualt: no filter
                                                if nargin < 4
                                                    true_thresh = [];
                                                    if nargin < 3
                                                        threshold = [0 1 2 0]; % default: 2 STD's from from absolute velocity median
                                                    end
                                                end
                                             end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            if struct_flag
                fnames = fieldnames(struct_inputs);
                for n = 1:length(fnames)
                    obj.(fnames{n}) = struct_inputs.(fnames{n});
                end
            else
                obj.sacd_length     = sacd_length;              % saccade length in time units
                obj.amp_cut         = amp_cut;                	% saccade amplitude cutoff
                obj.dur_cut         = dur_cut;               	% saccade duration cutoff
                obj.Fc_detect       = Fc_detect;               	% detection filter frequencies
                obj.Fc_ss           = Fc_ss;                 	% start-stop filter frequencies
                obj.threshold       = threshold;                % detetction threshold
                obj.true_thresh     = true_thresh;           	% cut-off for peak velocity in raw signal
                obj.direction     	= sign(direction);       	% saccade direction to detect (-1,0,1)
                obj.pks           	= pks;                      % saccade peak locations
                obj.min_pkdist     	= min_pkdist;           	% minimum peak distance
                obj.min_pkwidth   	= min_pkwidth;            	% minimum peak width
                obj.min_pkprom     	= min_pkprom;            	% minimum peak promiance
                obj.min_pkthresh  	= min_pkthresh;           	% minimum height difference between peak and its neighbors
                obj.boundThresh    	= boundThresh;           	% ratio of peaks to determine start and stop
            end
            
            obj.time          	= time(:);                	% time vector
            obj.n             	= length(obj.time);       	% # of data points
            obj.index       	= (1:obj.n)';             	% index vector
            obj.Ts           	= mean(diff(obj.time));    	% sampling time
            obj.Fs           	= 1/obj.Ts;               	% sampling frequency
            
            if all(size(position) > 1) % position is multidimensional
             	[obj.position, position_new_filt] = two2one_position(obj, position, [8 nan]); % convert
            else
                obj.position = position(:); % position vector
                position_new_filt = obj.position;
            end
            obj.velocity      	= central_diff(obj.position, obj.Ts); % velocity
          	obj.acceleration 	= central_diff(obj.velocity, obj.Ts); % acceleration
            obj.abs_velocity	= abs(obj.velocity); % absolute velocity
            
            % Filter position & velocity for peak detection
            [obj.position_filt_detect,obj.velocity_filt_detect ] = ....
                filter_lo_hi(obj, position_new_filt, obj.Fc_detect);
            
          	% Filter position & velocity for start-stop
            [obj.position_filt_ss, obj.velocity_filt_ss] = ...
                filter_lo_hi(obj, position_new_filt, obj.Fc_ss);
            
            % Calculate saccade detcetion threshold
            obj = calculateThreshold(obj, obj.threshold);

            % Detect, isolate, & normalize sacades & inter-saccade intervals
            obj = detectSaccades(obj, obj.pks, obj.min_pkdist, obj.min_pkwidth, ...
                                    obj.min_pkprom, obj.min_pkthresh, obj.boundThresh);
            obj = removeSaccades(obj);
            obj = normSaccade(obj);
            
            if showplot
                plotSaccade(obj)
                if obj.count~=0
                    plotInterval(obj)
                else
                    warning('No saccades detetected')
                end
            end
        end
        
    	function obj = setDeafults(obj)
            % setDeafults: sets values to default
            % 
            
            empty_fields = ["time","position","velocity"];
            
            for f = 1:length(empty_fields)
                obj.normpeak_saccade.(empty_fields(f))      = [];
                obj.norm_interval.(empty_fields(f))         = [];
                obj.normstart_interval.(empty_fields(f))   	= [];
                obj.normend_interval.(empty_fields(f))      = [];
                obj.normstart_stimulus.(empty_fields(f))  	= [];
                obj.normend_stimulus.(empty_fields(f))   	= [];
                obj.error.(empty_fields(f))                 = [];
                obj.int_error.(empty_fields(f))             = [];
            end
        end
        
        function [position_filt,velocity_filt] = filter_lo_hi(obj,position,Fc)
            % filter_lo_hi: low- & high-pass filter position and velocity data
            %   Fc   : 2x1 [lp, hp] cuttoff frequencies
            %
            
            position_filt = [];
            
            % Low-pass
            if ~isnan(Fc(1))
                [b_lo,a_lo] = butter(3, Fc(1) / (obj.Fs/2), 'low');
                position_filt = filtfilt(b_lo, a_lo, position);
                velocity_filt = central_diff(position_filt, obj.Ts);
                velocity_filt = filtfilt(b_lo, a_lo, velocity_filt);
            end
            
            % High-pass
          	if ~isnan(Fc(2))
                if isempty(position_filt) % no lp applied
                    temp_pos = position;
                else
                    temp_pos = position_filt;
                end
                [b_hi,a_hi] = butter(3, Fc(2) / (obj.Fs/2), 'high');
                position_filt = filtfilt(b_hi, a_hi, temp_pos);
              	velocity_filt = central_diff(position_filt, obj.Ts);
                velocity_filt = filtfilt(b_hi, a_hi, velocity_filt);
            end
            
            % No lp or hp filter applied
            if isempty(position_filt)
                position_filt = position;
              	velocity_filt = central_diff(position_filt, obj.Ts);
            end
        end
        
        function obj = calculateThreshold(obj, thresh)
            % calculateThreshold: compute detection threshold
            %   thresh   : [manual_thresh , med_boolean, nstd_thresh, min_thresh]
            %
            
            assert( length(thresh) <= 4, 'threshold must be 1x1 - 4x1')
            nT = length(thresh);
            
            obj.medians.velocity        = median(obj.velocity_filt_detect); % median velocity
            obj.stds.velocity       	= std(obj.velocity_filt_detect); % STD velocity
            obj.medians.abs_velocity	= median(abs(obj.velocity_filt_detect )); % median absolute velocity
            obj.stds.abs_velocity     	= std(abs((obj.velocity_filt_detect ))); % STD absolute velocity
            
            if nT >= 1
                thresh_all(1) = thresh(1);
                thresh_all(2) = -thresh(1);
                
                if nT >= 3 % calculate threshold based on # of standard deviations
                    assert(thresh(3) >= 0, 'STD threshold must be positive');
                    obj.nstd_thresh = thresh(3);
                    std_thresh = obj.nstd_thresh*obj.stds.velocity;
                else
                    std_thresh = 0;
                    obj.nstd_thresh = nan;
                end
                
                if nT >= 2
                    assert( any(thresh(2) == [0 1 2]), '2nd element of threshold must be 0, 1, or 2')
                    if thresh(2) == 0 % absolute threshold
                        temp_med = 0;
                        temp_thresh = thresh(1);
                        thresh_all(1) = temp_med + temp_thresh + std_thresh;
                        thresh_all(2) = temp_med - temp_thresh - std_thresh;
                    elseif thresh(2) == 1 % threshold from absolute median
                        temp_med = obj.medians.abs_velocity;
                        temp_thresh = thresh(1);
                        thresh_all(1) = temp_med + temp_thresh + std_thresh;
                        thresh_all(2) = -thresh_all(1);
                    elseif thresh(2) == 2 % threshold from raw median (only uneven threshold)
                        temp_med = obj.medians.velocity;
                        temp_thresh = thresh(1);
                        thresh_all(1) = temp_med + temp_thresh + std_thresh;
                        thresh_all(2) = temp_med - temp_thresh - std_thresh;
                    end
                end
            end
            
         	% Ensure minimum threshold
            if nT == 4
                if thresh_all(1) < thresh(4)
                   thresh_all(1) = thresh(4);
                end
                if thresh_all(2) > -thresh(4)
                   thresh_all(2) = -thresh(4);
                end
            end
            obj.threshold = thresh_all;
        end
        
        function [position_new,position_new_filt] = two2one_position(obj, position, Fc)
            % two2one_position: transforms nx2 position signal x to: y = x(:,1) - x(:,2)
            %   Fc   : optional cutoff frequency fro filtering
            %
            [~,dI] = max(size(position));
            assert(dI <= 2, 'position must be nx2 or 2xn')
            if dI == 2
                position = position';
            end
            position_new = position(:,1) - position(:,2);
            
            if nargin >= 2
                [position_filt_1,~] = filter_lo_hi(obj, position(:,1), Fc);
                [position_filt_2,~] = filter_lo_hi(obj, position(:,2), Fc);
                position_new_filt = position_filt_1 - position_filt_2;
            else
                position_new_filt = position;
            end
        end
        
        function obj = detectSaccades(obj, pks, min_pkdist, min_pkwidth, min_pkprom, min_pkthresh, boundThresh)
            % detetcSaccades: detetc saccades in data
            %   pks          : saccade peak locations [indicies] (optional)
            %   min_pkdist   : minimum peak distance [s]
            %   min_pkwidth  : minimum peak width [s]
            %   min_pkprom   : minimum peak promiance
            %   min_pkthresh : minimum height difference between peak and its neighbors
            %   boundThresh  : ratio of peaks to determine start and stop
            %
            
            % Set defaults
            if nargin < 7
                boundThresh = 0.25; % default bounding threshold
                if nargin < 6
                    min_pkthresh = 0; % defalt minimum height seperation threshold
                    if nargin < 5
                        min_pkprom = 0; % default peak prominence
                        if nargin < 4
                            min_pkwidth = 0; % default minimum peak width [s]
                            if nargin < 3
                                min_pkdist = 0; % default minimum peak distance [s]
                            end
                        end
                    end
                end
            end
            
            % Store detection properties
         	obj.boundThresh     = boundThresh;
            obj.min_pkthresh    = min_pkthresh;
            obj.min_pkprom      = min_pkprom;
            obj.min_pkwidth     = min_pkwidth;
            obj.min_pkdist      = min_pkdist;
            
            if isempty(pks) % if peaks are not given
                obj.ismanual = false; % automatically detetc saccades
                svel = obj.velocity_filt_detect; % detect saccade in this signal
                
                % Find peaks in filtered velocity data in both directions
                [~,locs_pos,~,p_pos] = findpeaks(svel, ...
                                        'MinPeakProminence', obj.min_pkprom, ...
                                        'MinPeakHeight', obj.threshold(1), ...
                                        'MinPeakDistance', round(obj.min_pkdist*obj.Fs), ...
                                        'MinPeakWidth', round(obj.min_pkwidth*obj.Fs), ...
                                        'Threshold', obj.min_pkthresh);
                [~,locs_neg,~,p_neg] = findpeaks(-svel, ...
                                        'MinPeakProminence', obj.min_pkprom, ...
                                        'MinPeakHeight', -obj.threshold(2), ...
                                        'MinPeakDistance', round(obj.min_pkdist*obj.Fs), ...
                                        'MinPeakWidth', round(obj.min_pkwidth*obj.Fs), ...
                                        'Threshold', obj.min_pkthresh);
                                    
             	% All peaks
               	[locs,locI] = sort([locs_pos ; locs_neg]);
                abs_peak = abs(svel(locs));
                %abs_peak = abs_peak(locI);
                dist = diff(locs); dist = [obj.n  ;dist];
                
                % Remove peaks to close together by picking largest peak
                %locs(dist < round(obj.min_pkdist*obj.Fs)) = [];
                locs_rmv = [];
                pp = 1;
                for ww = 1:length(locs)
                   if dist(ww) < round(obj.min_pkdist*obj.Fs) % too close
                       I_set = ww-1:ww;
                       peak_set = abs_peak(I_set);
                       [~,keepI] = min(peak_set);
                       locs_rmv(pp) = ww - length(I_set) + keepI;
                       pp =  pp + 1;
                   end
                end
                locs(locs_rmv) = [];

                % Exclude saccades in very beginning and end of signal
                window = round(obj.Fs*0.02); % window length at start & end to ignore saccades [samples]
                %window = 0;
                I = (locs > window) & (locs < (obj.n - (window+1))); % saccades inside middle window
                obj.peaks.index = locs(I);             
                
            else % user specified peaks
             	obj.ismanual = true; % don't auotmatically detect saccades
                if size(pks,2) == 3 % starts and ends are specified
                    obj.peaks.index = pks(:,2);
                elseif size(pks,2) == 1  % just peaks
                    obj.peaks.index = pks;
                else
                    error('Passed "peaks" must be nx1 or nx3')
                end
                    
                obj.threshold = [nan nan]; % no threshold for manual mode
                obj.Fc_detect = nan(1,2); % no detection filter needed
                %obj.velocity_filt_detect = nan(obj.n,1); % no detection filtered signal needed
            end
                     
            % Only use saccades in the direction specified
            if obj.direction==0
                % analyze all saccades
            else % remove saccades that don't match the input direction
                dir_include = sign(obj.velocity_filt_detect(obj.peaks.index)) == obj.direction;
                obj.peaks.index = obj.peaks.index(dir_include);
            end
            obj.count = length(obj.peaks.index); % # of saccades
            
            % Find true peaks in unfiltered data
            obj.peaks.index_detect = obj.peaks.index;
            peak_win = round(0.02*obj.Fs); % window to find true peaks
            for ww = 1:obj.count % each saccade peak
                Echeck = (obj.peaks.index(ww) - peak_win) : (obj.peaks.index(ww) + peak_win);
                Echeck = Echeck( (Echeck > 0) & (Echeck <= obj.n) );
                search_win = nan(obj.n,1);
                search_win(Echeck) = obj.velocity(Echeck);
                if obj.velocity_filt_detect(obj.peaks.index(ww)) >= 0
                        [~,shiftI] = max(search_win);
                elseif obj.velocity_filt_detect(obj.peaks.index(ww)) < 0
                        [~,shiftI]  = max(-search_win);
                end
                obj.peaks.index(ww) = shiftI;
            end
            
            % Don't include saccades if peak velocity is less than true_thresh
            if ~isempty(obj.true_thresh)
                above_thresh = abs(obj.velocity(obj.peaks.index)) >= obj.true_thresh;
                obj.peaks.index = obj.peaks.index(above_thresh);
            end
            
            % Get peak properties for detecting start/end of saccades            
            obj.peaks.velocity	= obj.velocity(obj.peaks.index); % velocity peaks
            obj.count          	= length(obj.peaks.index); % # of saccades
            
            % Only for detection
          	obj.peaks.velocity_detect   = obj.velocity_filt_detect(obj.peaks.index); % detection velocity peaks
            obj.peaks.velocity_ss       = obj.velocity_filt_ss(obj.peaks.index); % start/stop velocity peaks
            
            % Names for saccade statistics table
            obj.tablenames = {  'Amplitude'  , 'Direction'   , 'Duration'   , ...
                                'RiseTime'  ,   'FallTime'   , 'Skew'       , ...
                                'StartIdx'  ,   'PeakIdx'    , 'EndIdx'     , ...
                                'StartTime' ,   'PeakTime'   , 'EndTime'    , ...
                                'StartPos'  ,   'PeakPos'    , 'EndPos'     , ...
                                'StartVel'  ,   'PeakVel'    , 'EndVel'};
              
            if ~isempty(obj.peaks.index) % if any saccades are detected   
                if isnan(obj.sacd_length) % automatically find start and ends of saccades
                    if size(pks,2) == 3  % starts and ends are specified
                        obj.starts.index = pks(:,1);
                        obj.peaks.index = pks(:,2);
                        obj.ends.index = pks(:,3);
                    else
                        for ww = 1:obj.count % every saccade
                            dir_sign = sign(obj.peaks.velocity(ww));
                            dir_positive = obj.peaks.velocity(ww) >= 0; % is the saccade positive or negative
                            bound = obj.peaks.velocity_ss(ww)*obj.boundThresh(1);
                            if length(obj.boundThresh) == 2
                                if abs(bound) < obj.boundThresh(2)
                                    if bound == 0
                                        sgn = 1;
                                    else
                                        sgn = sign(bound);
                                    end
                                    bound = sgn * obj.boundThresh(2);
                                end
                            end
                            
                            % Find start of saccade
                            if dir_positive
                                Sidx = find(obj.velocity_filt_ss(1:obj.peaks.index(ww)) ...
                                                        <= bound, 1,'last');
                            else
                                Sidx = find(obj.velocity_filt_ss(1:obj.peaks.index(ww)) ...
                                                        >= bound, 1,'last');
                            end

                            % Make sure we have the start of the saccade
                            if ~isempty(Sidx)
                                % Make sure saccade starts with velocity away from 0
                                check_span = [Sidx , Sidx+1];
                                check_dir = sign(diff(obj.velocity(check_span))) ~= dir_sign;
                                if check_dir
                                    Scheck = 1:Sidx;
                                    Ss = find(dir_sign*obj.velocity(Scheck) < ...
                                        dir_sign*obj.velocity(Sidx), 1, 'last');
                                    if ~isempty(Ss)
                                        Sidx = Scheck(Ss);
                                    end
                                end
                                
                                obj.starts.index(ww,1) = Sidx;
                            else
                                obj.starts.index(ww,1) = nan;
                            end

                            % Find end of saccade
                            if dir_positive
                                Echeck = find(obj.velocity_filt_ss <= bound);
                            else
                                Echeck = find(obj.velocity_filt_ss >= bound);
                            end
                            Es = find(Echeck > obj.peaks.index(ww),1,'first');
                            Eidx = Echeck(Es);
                                                        
                            % Make sure we have the end of the saccade
                            if ~isempty(Es) % make sure saccade ends
                                % Make sure saccade ends with velocity towards 0
                                check_span = [Eidx-1 , Eidx];
                                check_dir = sign(diff(obj.velocity(check_span))) == dir_sign;
                                if check_dir
                                    Echeck = Eidx-1:obj.n;
                                    Es = find(dir_sign*obj.velocity(Echeck) < ...
                                        dir_sign*obj.velocity(Eidx), 1, 'first');
                                    if ~isempty(Es)
                                        Eidx = Echeck(Es);
                                    end
                                end
                                
                                obj.ends.index(ww,1) = Eidx; % saccade end index
                            else
                                obj.ends.index(ww,1) = nan; % saccade end index not found
                            end
                        end
                    end
                else % start and ends of saccades are set distance from the peak time
                    Echeck = ceil((obj.sacd_length * obj.Fs)/2);
                    obj.starts.index = obj.peaks.index - Echeck;
                    obj.ends.index = obj.peaks.index + Echeck;
                    obj.ends.index(obj.ends.index > obj.n) = obj.n;
                    obj.starts.index(obj.starts.index < 1) = 1;
                end
                
                % If saccade is missing a start or end, don't include it
                missing_start_end   = logical(sum(isnan([obj.starts.index, obj.ends.index]),2));
                obj.starts.index    = obj.starts.index(~missing_start_end);
                obj.peaks.index     = obj.peaks.index(~missing_start_end);
                obj.ends.index      = obj.ends.index(~missing_start_end);
                obj.peaks.velocity  = obj.peaks.velocity(~missing_start_end);
                obj.count           = length(obj.peaks.index);

                % Get peak, start, and end values
                obj.peaks.time   	= obj.time(obj.peaks.index); % time peaks
                obj.peaks.position  = obj.position(obj.peaks.index); % position peaks
                
                obj.starts.time   	= obj.time      (obj.starts.index);
                obj.starts.position	= obj.position  (obj.starts.index);
                obj.starts.velocity	= obj.velocity  (obj.starts.index);
                
                obj.ends.time       = obj.time      (obj.ends.index);
                obj.ends.position   = obj.position  (obj.ends.index);
                obj.ends.velocity   = obj.velocity  (obj.ends.index);
                
                % Saccade amplitudes % durations
                obj.amplitudes = obj.ends.position - obj.starts.position;
                obj.durations  = obj.ends.time - obj.starts.time;
                
                % Remove saccades below position amplitude threshold & duration threshold
              	rmv_amp = abs(obj.amplitudes) < obj.amp_cut;
                rmv_dur = obj.durations > obj.dur_cut;
                rmv_all = rmv_amp | rmv_dur;
                if any(rmv_all)
                    obj.amplitudes = obj.amplitudes(~rmv_all,:);

                    obj.peaks.index  	= obj.peaks.index(~rmv_all);
                    obj.peaks.time      = obj.peaks.time(~rmv_all);
                    obj.peaks.position  = obj.peaks.position(~rmv_all);
                    obj.peaks.velocity  = obj.peaks.velocity(~rmv_all);
                    
                    obj.starts.index  	= obj.starts.index(~rmv_all);
                    obj.starts.time   	= obj.starts.time(~rmv_all);
                    obj.starts.position = obj.starts.position(~rmv_all);
                    obj.starts.velocity = obj.starts.velocity(~rmv_all);

                    obj.ends.index  	= obj.ends.index(~rmv_all);
                    obj.ends.time       = obj.ends.time(~rmv_all);
                    obj.ends.position   = obj.ends.position(~rmv_all);
                    obj.ends.velocity   = obj.ends.velocity(~rmv_all);
                    
                    obj.durations     	= obj.ends.time - obj.starts.time;

                    obj.count           = length(obj.peaks.index);
                end

                rise_time           = obj.peaks.time - obj.starts.time;
                fall_time           = obj.ends.time - obj.peaks.time;
                skew                = rise_time ./ obj.durations;
                obj.directions      = sign(obj.peaks.velocity);
             	obj.rate            = obj.count / (obj.time(end) - obj.time(1));
                obj.cmap            = hsv(obj.count);
                
                obj.peaks.velocity_detect   = obj.velocity_filt_detect(obj.peaks.index); % detection velocity peaks
                obj.peaks.velocity_ss       = obj.velocity_filt_ss(obj.peaks.index); % start/stop velocity peaks
                
                if obj.count == 0
                    STATS = nan(size(obj.tablenames));
                else
                    % Collect data
                    STATS = [obj.amplitudes       , obj.directions    	, obj.durations         , ...
                             rise_time            , fall_time           , skew                  , ...
                             obj.starts.index     , obj.peaks.index     , obj.ends.index        , ...
                             obj.starts.time	  , obj.peaks.time      , obj.ends.time         , ...
                             obj.starts.position  , obj.peaks.position  , obj.ends.position     , ...
                             obj.starts.velocity  , obj.peaks.velocity  , obj.ends.velocity     ];
                end
            	
                % Saccade table
                obj.SACD = splitvars(table(STATS));
                obj.SACD.Properties.VariableNames = obj.tablenames;
                
            else
                STATS = nan(size(obj.tablenames));
                obj.count = 0;
                obj.rate = 0;
                
             	% Saccade table
                obj.SACD = splitvars(table(STATS));
                obj.SACD.Properties.VariableNames = obj.tablenames;
            end
            
        end
        
      	function obj = removeSaccades(obj)
      	% removeSaccades: extract saccade kinematics from data & grab
      	% intervals between saccades
     	% 
            if obj.count~=0 % if any saccades are detected
                obj.saccades  = cell(obj.count,1); % pull out each saccade
                obj.intervals = cell(obj.count,1); % pull out each interval
                for ww = 1:obj.count % every saccade
                    % Create saccade table
                    n_saccade_idx = obj.SACD.EndIdx(ww) - obj.SACD.StartIdx(ww) + 1;
                    obj.saccades{ww} = splitvars(table(nan(n_saccade_idx,4)));
                    obj.saccades{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                    
                    % Assign variables
                    obj.saccades{ww}.Index      = (obj.SACD.StartIdx(ww):obj.SACD.EndIdx(ww))';
                 	obj.saccades{ww}.Time   	= obj.time(obj.saccades{ww}.Index);
                    obj.saccades{ww}.Position 	= obj.position(obj.saccades{ww}.Index);
                    obj.saccades{ww}.Velocity 	= obj.velocity(obj.saccades{ww}.Index);
                    
                    % Create interval table
                    if ww~=1
                        n_interval_idx = obj.SACD.StartIdx(ww) - obj.SACD.EndIdx(ww-1) + 1;
                        if n_interval_idx < 1
                            obj.intervals{ww} = splitvars(table(nan(3,4)));
                            obj.intervals{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                        else
                            obj.intervals{ww} = splitvars(table(nan(n_interval_idx,4)));
                            obj.intervals{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                            obj.intervals{ww}.Index     = (obj.SACD.EndIdx(ww-1):obj.SACD.StartIdx(ww))';
                            obj.intervals{ww}.Time      = obj.time(obj.intervals{ww}.Index);
                            obj.intervals{ww}.Position  = obj.position(obj.intervals{ww}.Index);
                            obj.intervals{ww}.Velocity  = obj.velocity(obj.intervals{ww}.Index);
                        end
                    else % store 1st interval as nan's because we don't know the start of this interval
                        n_interval_idx = obj.SACD.StartIdx(ww);
                        if n_interval_idx < 3
                           n_interval_idx = 3; 
                        end
                        obj.intervals{ww} = splitvars(table(nan(n_interval_idx,4)));
                     	obj.intervals{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                        
                        obj.intervals{ww}.Index     = (1:n_interval_idx)';
                        obj.intervals{ww}.Time      = nan*obj.time(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Position  = nan*obj.position(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Velocity  = nan*obj.velocity(obj.intervals{ww}.Index);
                    end
                end
                
                % Get only saccade data
                obj.saccades_all = cat(1,obj.saccades{:});
                saccade_idx = ismember(obj.index,obj.saccades_all.Index)';
                obj.sacd_ratio = sum(saccade_idx) / obj.n;
                
                % Remove saccades & replace with nan's
                obj.removed_all = table(obj.index,obj.time,obj.position,obj.velocity,...
                    'VariableNames',{'Index','Time','Position','Velocity'});
                
                obj.removed_all.Index(saccade_idx)      = nan;
                obj.removed_all.Time(saccade_idx)       = nan;
                obj.removed_all.Position(saccade_idx)   = nan;
                obj.removed_all.Velocity(saccade_idx)   = nan;
                
                % Setup table to store shifted data with removed saccades
                obj.shift = obj.removed_all;
                obj.shift.Index = obj.index;
                obj.shift.Time = obj.time;
                
                % Shift position between saccades so data after a saccade
                % starts where the last saccades ends
                for ww = 1:obj.count % every saccade
                    obj.shift.Position(obj.SACD.EndIdx(ww):end) = obj.shift.Position(obj.SACD.EndIdx(ww):end) ...
                        - ( obj.SACD.EndPos(ww) - obj.SACD.StartPos(ww) );
                end
                
                % Remove nan's for interpolation
                time_shift = obj.removed_all.Time(~isnan(obj.removed_all.Time));
                pos_shift  = obj.shift.Position(~isnan(obj.removed_all.Time));
                %vel_shift  = obj.shift.Velocity(~isnan(obj.removed_all.Time));
                
                % Interpolate between saccades
                intrp_pos = table( interp1(time_shift, pos_shift, obj.shift.Time, 'spline'),... % interpolate
                    'VariableNames',{'IntrpPosition'} );
                intrp_vel = central_diff(intrp_pos{:,1}, obj.Ts); % velocity
                intrp_vel = table( intrp_vel, 'VariableNames',{'IntrpVelocity'} );
            else
                intrp_pos = table( obj.position,'VariableNames',{'IntrpPosition'} );
                intrp_vel = table( obj.velocity,'VariableNames',{'IntrpVelocity'} );
            end
            % Add to shift table
            obj.shift = [obj.shift, intrp_pos, intrp_vel];
        end
        
        function [obj] = normSaccade(obj)
        % normSaccade: align saccades to peak time & align intervals to
        % start and end time.
        %  
        
            int_varnames = {'IntTime','IntAmp','IntRange','IntMeanVel'};
            int_after_varnames = cellstr(string(int_varnames) + "_after");
            if obj.count~=0 % if there are any saccades
                % Normalize saccade times to saccade peak times
                obj.normpeak_saccade.time = cellfun(@(x,y) x.Time - y, obj.saccades, ...
                    num2cell(obj.SACD.PeakTime),'UniformOutput', false);

                [obj.normpeak_saccade.time ,~,~,~,obj.align_saccade_peak,~] = ...
                    nancat_center(obj.normpeak_saccade.time, 0, 1);

                % Align saccade positions to saccade peak times
                obj.normpeak_saccade.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), obj.saccades, ...
                    obj.align_saccade_peak, 'UniformOutput', false);
                obj.normpeak_saccade.position = cat(2,obj.normpeak_saccade.position{:});

                % Align saccade velocities to saccade peak times
                obj.normpeak_saccade.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.saccades, ...
                    obj.align_saccade_peak, 'UniformOutput', false);
                obj.normpeak_saccade.velocity = cat(2,obj.normpeak_saccade.velocity{:});

                % Compute saccade stats
                obj.normpeak_saccade.time_stats = basic_stats(obj.normpeak_saccade.time,2);
                obj.normpeak_saccade.position_stats = basic_stats(obj.normpeak_saccade.position,2);
                obj.normpeak_saccade.velocity_stats = basic_stats(obj.normpeak_saccade.velocity,2);

                %-------------------------------------------

                % Normalize interval times to intervals start times
                obj.normstart_interval.time = cellfun(@(x) x.Time - x.Time(1), obj.intervals, ...
                    'UniformOutput', false);

                [obj.normstart_interval.time ,~,~,~,obj.align_interval_start,~] = ...
                    nancat_center(obj.normstart_interval.time, 0 ,1);
               
                % Align interval start times without normalizing positions
                obj.norm_interval.time = obj.normstart_interval.time;
                
                obj.norm_interval.position = cellfun(@(x,y) padmat(x.Position, y,nan,1), ...
                    obj.intervals, obj.align_interval_start, 'UniformOutput', false);
                obj.norm_interval.position = cat(2,obj.norm_interval.position{:});

                % Align interval positions to interval start times
                obj.normstart_interval.position = cellfun(@(x,y) padmat(x.Position - x.Position(1),y,nan,1), ...
                    obj.intervals, obj.align_interval_start, 'UniformOutput', false);
                obj.normstart_interval.position = cat(2,obj.normstart_interval.position{:});

                % Align interval velocities to interval start times
                obj.normstart_interval.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.intervals, ...
                    obj.align_interval_start, 'UniformOutput', false);
                obj.normstart_interval.velocity = cat(2,obj.normstart_interval.velocity{:});
                
                obj.norm_interval.velocity = obj.normstart_interval.velocity;

                % Compute interval stats
                obj.normstart_interval.time_stats = basic_stats(obj.normstart_interval.time,2);
                obj.normstart_interval.position_stats = basic_stats(obj.normstart_interval.position,2);
                obj.normstart_interval.velocity_stats = basic_stats(obj.normstart_interval.velocity,2);
                
                obj.norm_interval.position_stats = basic_stats(obj.norm_interval.position,2);

                % Get interval times, ampltitudes, & range
                obj.normstart_interval.endidx = sum(~isnan(obj.normstart_interval.time));
                
                zeroIdx = obj.normstart_interval.endidx(1,:) == 0;
                obj.normstart_interval.endidx(zeroIdx) = obj.normstart_interval.endidx(zeroIdx) + 1;
                
                IntTime = max(obj.normstart_interval.time(obj.normstart_interval.endidx,:))';
                IntRange = range(obj.norm_interval.position,1)';
                IntMeanVel = nanmean(obj.norm_interval.velocity,1)';
                
                IntAmp = nan(obj.count,1);
                for ww = 1:obj.count
                    IntAmp(ww,1) = obj.normstart_interval.position(obj.normstart_interval.endidx(ww),ww);
                end
                
                Int_All = table(IntTime, IntAmp, IntRange, IntMeanVel, ...
                    'VariableNames', int_varnames);
                
                Int_All_after = circshift(Int_All, -1, 1);
                Int_All_after.Properties.VariableNames = int_after_varnames;
                
                obj.SACD = [obj.SACD , Int_All, Int_All_after];

                %-------------------------------------------

                % Normalize interval times to intervals end times
                obj.normend_interval.time = cellfun(@(x) x.Time - x.Time(end), obj.intervals, ...
                    'UniformOutput', false);

                [obj.normend_interval.time ,~,~,~,obj.align_interval_end,~] = ...
                    nancat_center(obj.normend_interval.time, 0 ,1);

                % Align interval positions to interval end times
                obj.normend_interval.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), obj.intervals, ...
                    obj.align_interval_end, 'UniformOutput', false);
                obj.normend_interval.position = cat(2,obj.normend_interval.position{:});

                % Align interval velocities to interval end times
                obj.normend_interval.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.intervals, ...
                    obj.align_interval_end, 'UniformOutput', false);
                obj.normend_interval.velocity = cat(2,obj.normend_interval.velocity{:});

                % Compute interval stats
                obj.normend_interval.time_stats = basic_stats(obj.normstart_interval.time,2);
                obj.normend_interval.position_stats = basic_stats(obj.normstart_interval.position,2);
                obj.normend_interval.velocity_stats = basic_stats(obj.normstart_interval.velocity,2);
            else
                all_varnames = [int_varnames, int_after_varnames];
                Int_All = splitvars(table(nan(size(all_varnames))));
             	Int_All.Properties.VariableNames = all_varnames;
                obj.SACD = [obj.SACD , Int_All];
                obj = setDeafults(obj);
            end
        end
        
        function [obj] = stimSaccade(obj,stim,debug)
        % stimSaccade: computes saccade/interval relationship to visual
        % motion stimulus (stimulus position, stimulus velocity, error,integrated error)
        %  
        
         	% Get stimulus data
            obj.stimlus_position = stim(:);
            obj.stimlus_velocity = central_diff(obj.stimlus_position, obj.Ts);
            obj.stimulus = cell(obj.count,1); % pull out each saccade
        
            if obj.count~=0 % if there are any saccades
                if nargin < 3
                    debug = false; % default
                end

                % Extract stimulus intervals
                vel_window = nan(obj.count,1); % stim velocity in window before saccade
                med_win_co_anti = round(0.5 * obj.Fs); % mean stimulus window size
                if ~isempty(stim) % if any saccades are detected and a stimulus if given
                    for ww = 1:obj.count % every saccade
                        % Create stimulus table
                        if ww~=1
                            n_stim_idx = obj.SACD.StartIdx(ww) - obj.SACD.EndIdx(ww-1) + 1;
                            if n_stim_idx > 0
                                obj.stimulus{ww} = splitvars(table(nan(n_stim_idx,4)));
                                obj.stimulus{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                                obj.stimulus{ww}.Index     = (obj.SACD.EndIdx(ww-1):obj.SACD.StartIdx(ww))';
                                obj.stimulus{ww}.Time      = obj.time(obj.stimulus{ww}.Index);
                                obj.stimulus{ww}.Position  = obj.stimlus_position(obj.stimulus{ww}.Index);
                                obj.stimulus{ww}.Velocity  = obj.stimlus_velocity(obj.stimulus{ww}.Index);
                            else
                                warning('Saccades too close together to find intervals')
                                obj.stimulus{ww} = splitvars(table(nan(3,4)));
                                obj.stimulus{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                                obj.stimulus{ww}.Index     = (1:4)';
                                obj.stimulus{ww}.Time      = nan*obj.time(obj.stimulus{ww}.Index);
                                obj.stimulus{ww}.Position  = nan*obj.stimlus_position(obj.stimulus{ww}.Index);
                                obj.stimulus{ww}.Velocity  = nan*obj.stimlus_velocity(obj.stimulus{ww}.Index);
                            end

                        else % store 1st interval as nan's because we don't know the start of this interval
                            n_stim_idx = obj.SACD.StartIdx(ww);
                            if n_stim_idx < 3
                                n_stim_idx = 3;
                            end
                            
                            obj.stimulus{ww} = splitvars(table(nan(n_stim_idx,4)));
                            obj.stimulus{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                            obj.stimulus{ww}.Index     = (1:n_stim_idx)';
                            obj.stimulus{ww}.Time      = nan*obj.time(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Position  = nan*obj.stimlus_position(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Velocity  = nan*obj.stimlus_velocity(obj.stimulus{ww}.Index);
                        end
                        vel_temp = obj.stimlus_velocity(obj.stimulus{ww}.Index);
                        int_length = length(obj.stimulus{ww}.Velocity);
                        if int_length >= med_win_co_anti
                            vel_window(ww) = mean(vel_temp(1 + int_length-med_win_co_anti:int_length));
                        else
                            vel_window(ww) = mean(obj.stimulus{ww}.Velocity);
                        end
                    end
                    
                    % Classify saccades as anti- or co-directional based on stimulus velocity
                    co_anti = obj.SACD.Direction .* sign(vel_window);
                    
                    % Assign times
                    obj.normstart_stimulus.time = obj.normstart_interval.time;
                    obj.normend_stimulus.time   = obj.normend_interval.time;
                    %-------------------------------------------
                    
                    % Align stimulus interval positions to interval start times
                    obj.normstart_stimulus.position = cellfun(@(x,y) padmat(x.Position - x.Position(1),y,nan,1), ...
                        obj.stimulus, obj.align_interval_start, 'UniformOutput', false);
                    obj.normstart_stimulus.position = cat(2,obj.normstart_stimulus.position{:});

                    % Align stimulus interval velocities to interval start times
                    obj.normstart_stimulus.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), ...
                        obj.stimulus, obj.align_interval_start, 'UniformOutput', false);
                    obj.normstart_stimulus.velocity = cat(2,obj.normstart_stimulus.velocity{:});
                    
                    %-------------------------------------------
                    
                    % Align stimulus interval positions to interval end times
                    obj.normend_stimulus.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), ...
                        obj.stimulus, obj.align_interval_end, 'UniformOutput', false);
                    obj.normend_stimulus.position = cat(2,obj.normend_stimulus.position{:});

                    % Align stimulus interval velocities to interval end times
                    obj.normend_stimulus.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), ...
                        obj.stimulus, obj.align_interval_end, 'UniformOutput', false);
                    obj.normend_stimulus.velocity = cat(2,obj.normend_stimulus.velocity{:});
                    
                    %-------------------------------------------

                    % Get time vectors for integration & compute error & integrated error within intervals
                    [~,max_time] = max(obj.normstart_interval.endidx);
                    tvector = obj.normstart_interval.time(:,max_time);

                    % Compute error within intervals
                    obj.error.time  	= obj.normstart_interval.time;
                    obj.error.position  = obj.normstart_stimulus.position - obj.normstart_interval.position;
                    obj.error.velocity  = obj.normstart_stimulus.velocity - obj.normstart_interval.velocity;

                    obj.error.position_end = arrayfun(@(x,y)  obj.error.position(y,x), ...
                                1:size(obj.error.position,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.error.velocity_end = arrayfun(@(x,y)  obj.error.velocity(y,x), ...
                                1:size(obj.error.velocity,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);

                    % Compute integrated error within intervals
                    obj.int_error.time      = obj.normstart_interval.time;
                    obj.int_error.position  = cumtrapz(tvector, obj.error.position, 1);
                    obj.int_error.velocity  = cumtrapz(tvector, obj.error.velocity, 1);

                    obj.int_error.position_end = arrayfun(@(x,y)  obj.int_error.position(y,x), ...
                                1:size(obj.int_error.position,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.int_error.velocity_end = arrayfun(@(x,y)  obj.int_error.velocity(y,x), ...
                                1:size(obj.int_error.velocity,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.int_error.position_end(1) = nan; % integrated error is nan for 1st saccade
                    obj.int_error.velocity_end(1) = nan; % integrated error is nan for 1st saccade

                    % Make table of finals erros & combine with SACD table
                    ERR = table(obj.error.position_end',     obj.error.velocity_end', ...
                                obj.int_error.position_end', obj.int_error.velocity_end', vel_window, co_anti, ...
                                'VariableNames', {'ErrorPos','ErrorVel','IntErrorPos','IntErrorVel',...
                                'StimMedVel', 'CoAnti'});
                    obj.SACD = [obj.SACD , ERR];

                    % Compute interval stats
                    obj.error.position_stats = basic_stats(obj.error.position,2);
                    obj.error.velocity_stats = basic_stats(obj.error.velocity,2);

                    obj.int_error.position_stats = basic_stats(obj.int_error.position,2);
                    obj.int_error.velocity_stats = basic_stats(obj.int_error.velocity,2);
                end

                % Plot if specified
                if debug
                    if ~obj.count==0
                        plotStimulus(obj) 
                    else
                        warning('No saccades detetected')
                    end
                end
                
            else
                % Saccade table
                err_table = table(nan, nan, nan, nan, nan, nan, 'VariableNames', ...
                     {'ErrorPos','ErrorVel','IntErrorPos','IntErrorVel','StimMedVel','CoAnti'});
                obj.SACD = [obj.SACD , err_table];
            end
        end
        
        function [scds,ints,scd_time,int_time] = getSaccade(obj, new_signal, scd_win, int_win, norm)
            % getSaccade: use start/stop times of saccades to pullout saccade and interval data from an 
            %             external signal
            %
            %   new_signal: new signal
            %   scd_win:    optional window around each saccade peak to pullout data. If not given then just use
            %               saccade start/end points
            %   int_win:    optional window from start of each interval peak to pullout data. If not given then just use
            %               saccade start/end points
            %   norm:       normalize saccades to be in the same direction
            %  
            
            if nargin < 5
                norm = [];
                if nargin < 4
                    int_win = [];
                    if nargin < 3
                        scd_win = [];
                    end
                end
            end
            
            if ~isempty(norm)
                dir = obj.SACD.Direction;
            else
                dir = ones(obj.count,1);
            end
            
            assert(obj.n == length(new_signal), 'Signal must be same length as saccade object''s signal')
            scd_winI = round(obj.Fs * scd_win); % saccade window size in samples
            int_winI = round(obj.Fs * int_win); % interval window size in samples
            scds = cell(obj.count,1);
            ints = cell(obj.count,1);
          	scd_time = cell(obj.count,1);
            int_time = cell(obj.count,1);
            for ww = 1:obj.count
                if isempty(scd_win)
                    scdI = obj.saccades{ww}.Index; % start to stop of each saccade
                    scds{ww} = new_signal(scdI);
                    scd_time{ww} = ( (scdI - obj.SACD.PeakIdx(ww)) * obj.Ts )';
                else
                    I = (obj.SACD.PeakIdx(ww) - scd_winI) : (obj.SACD.PeakIdx(ww) + scd_winI);
                   	outI = (I < 1) | (I > obj.n);
                    scdI = I(~outI);
                    scd_time{ww} = ( (I - obj.SACD.PeakIdx(ww)) * obj.Ts )';
                    scd_time{ww}(outI) = nan;
                    scds{ww} = nan(length(I),1);
                    scds{ww}(~outI) = dir(ww) * new_signal(scdI);
                end

                if ww == 1
                    nanflag = nan;
                else
                    nanflag = 1;
                end
                
                if isempty(int_win)
                    intI = obj.intervals{ww}.Index; % start to stop of each interval
                    ints{ww} = nanflag * new_signal(intI);
                    int_time{ww} = obj.Ts * (0:length(ints{ww})-1)';
                else
                    I = obj.intervals{ww}.Index(1) : (obj.intervals{ww}.Index(1) + int_winI);
                    if ~isnan(I)
                        outI = (I < 1) | (I > obj.n);
                        intI = I(~outI);
                        int_time{ww} = ( (I - obj.intervals{ww}.Index(1)) * obj.Ts )';
                        int_time{ww}(outI) = nan;
                        ints{ww} = nan(length(I),1);
                        ints{ww}(~outI) = dir(ww) * nanflag * new_signal(intI);
                    end
                end
            end
        end
        
        function plotSaccade(obj)
            % plotSaccade: plots data & detected saccades
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Saccades';
            set(FIG,'WindowStyle','docked')

            ax(1) = subplot(4,1,1) ; hold on
                ylabel('Position')
                h(1) = plot(obj.time,obj.position,'k');
                h_filt(1,1) = plot(obj.time,obj.position_filt_ss, 'Color' , [0.5 0.5 0.5]);
                plot(obj.time, zeros(obj.n,1),'--','Color',[0.5 0.5 0.5])
                if obj.count~=0
                    for ww = 1:obj.count
                       plot(obj.saccades{ww}.Time,obj.saccades{ww}.Position,...
                           'LineWidth', 1, 'Color', obj.cmap(ww,:))
                       plot(obj.intervals{ww}.Time,obj.intervals{ww}.Position,...
                           'LineWidth', 1, 'Color', 0.7*obj.cmap(ww,:))
                    end
                    plot(obj.starts.time , obj.starts.position , '*g')
                    plot(obj.peaks.time  , obj.peaks.position  , '*b')
                    plot(obj.ends.time   , obj.ends.position   , '*r')
                end
                h_filt(1,2) = plot(obj.time,obj.position_filt_detect,'k');
                %ax(1).YLim = max(abs(ax(1).YLim))*[-1 1];

            ax(2) = subplot(4,1,2) ; hold on
                ylabel('Velocity')
                h(2) = plot(obj.time,obj.velocity,'k');
                h_filt(2,1) = plot(obj.time,obj.velocity_filt_ss, 'Color' , [0.5 0.5 0.5]);
                plot(obj.time, obj.threshold(1)*ones(obj.n,1) , '--m')
                plot(obj.time,  obj.threshold(2)*ones(obj.n,1) , '--m')
                if obj.count~=0
                    for ww = 1:obj.count
                       plot(obj.saccades{ww}.Time,obj.saccades{ww}.Velocity,...
                           'LineWidth', 1, 'Color', obj.cmap(ww,:))
                       plot(obj.intervals{ww}.Time,obj.intervals{ww}.Velocity,...
                           'LineWidth', 1, 'Color', 0.7*obj.cmap(ww,:))
                    end
                    plot(obj.starts.time , obj.starts.velocity  , '*g')
                    plot(obj.peaks.time  , obj.peaks.velocity   , '*b')
                    plot(obj.ends.time   , obj.ends.velocity    , '*r')
                end
                h_filt(2,2) = plot(obj.time,obj.velocity_filt_detect,'k');
                ax(2).YLim = max(abs(ax(2).YLim))*[-1 1];

            ax(3) = subplot(4,1,3) ; hold on
                if obj.count~=0
                    ylabel('Removed Position')                
                    plot(obj.shift.Time, obj.shift.IntrpPosition,      	 'r', 'LineWidth', 1);
                    plot(obj.shift.Time, obj.shift.Position,             'k', 'LineWidth', 1);
                    plot(obj.removed_all.Time, obj.removed_all.Position, 'c', 'LineWidth', 1);
                end
            ax(4) = subplot(4,1,4) ; hold on
                if obj.count~=0
                    ylabel('Removed Velocity')
                    xlabel('Time')
                    plot(obj.shift.Time, obj.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
                    plot(obj.shift.Time, obj.shift.Velocity,        'k', 'LineWidth', 0.5);
                end
            set(ax,'LineWidth',1,'FontWeight','bold')
            linkaxes(ax,'x')
            set(h,'LineWidth',0.5)
            set(h_filt, 'LineWidth', 0.5)
            %align_Ylabels(FIG)
        end
        
        function plotInterval(obj)
            % plotInterval: plots extracted normalized saccades & intervals
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Normalized Saccades & Intervals';
            set(FIG,'WindowStyle','docked')
            ax(1) = subplot(2,2,1) ; hold on ; title('Saccades')
                ylabel('Position')
                h.sacdpos = plot(obj.normpeak_saccade.time,obj.normpeak_saccade.position);
                [hstd(1),~] = PlotPatch(obj.normpeak_saccade.position_stats.median, ...
                                obj.normpeak_saccade.position_stats.std, ...
                                obj.normpeak_saccade.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(1),'bottom')
                                
          	ax(2) = subplot(2,2,3) ; hold on
                ylabel('Velocity')
                h.sacdvel = plot(obj.normpeak_saccade.time,obj.normpeak_saccade.velocity);
               	[hstd(2),~] = PlotPatch(obj.normpeak_saccade.velocity_stats.median, ...
                                obj.normpeak_saccade.velocity_stats.std, ...
                                obj.normpeak_saccade.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(2),'bottom')
                                
          	ax(3) = subplot(2,2,2) ; hold on ; title('Intervals')
                h.intpos = plot(obj.normstart_interval.time,obj.norm_interval.position);
                [hstd(3),~] = PlotPatch(obj.norm_interval.position_stats.median, ...
                                obj.norm_interval.position_stats.std, ...
                                obj.normstart_interval.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(3),'bottom')
                            
          	ax(4) = subplot(2,2,4) ; hold on
                h.intvel = plot(obj.normstart_interval.time,obj.normstart_interval.velocity);
                [hstd(4),~] = PlotPatch(obj.normstart_interval.velocity_stats.median, ...
                                obj.normstart_interval.velocity_stats.std, ...
                                obj.normstart_interval.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(4),'bottom')
                            
            set(h.sacdpos, {'color'}, num2cell(obj.cmap,2))
            set(h.sacdvel, {'color'}, num2cell(obj.cmap,2))
            set(h.intpos,  {'color'}, num2cell(obj.cmap,2))
            set(h.intvel,  {'color'}, num2cell(obj.cmap,2))
            
          	set([h.sacdpos h.sacdvel h.intpos h.intvel],'LineWidth',1)
            set(ax,'LineWidth',1,'FontWeight','bold')
            linkaxes(ax(1:2),'x')
            linkaxes(ax(3:4),'x')
            %align_Ylabels(FIG)
        end
        
        function plotStimulus(obj)
            % plotStimulus: plots extracted stimulus, error, & integrated error intervals
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Stimulus & Error';
            set(FIG,'WindowStyle','docked')
            if obj.count > 1
                ax(1) = subplot(2,3,1) ; hold on ; title('Stimulus')
                    ylabel('Position')
                    h.stimpos = plot(obj.normstart_interval.time, obj.normstart_stimulus.position);

                ax(2) = subplot(2,3,4) ; hold on
                    ylabel('Velocity')
                    h.stimvel = plot(obj.normstart_interval.time, obj.normstart_stimulus.velocity);
                    ax(2).YLim = 1.1*abs(max(max(obj.normstart_stimulus.velocity)))*[-1 1];

                ax(3) = subplot(2,3,2) ; hold on ; title('Error')
                    h.errpos = plot(obj.normstart_interval.time,obj.error.position);
                    [hstd(1),~] = PlotPatch(obj.error.position_stats.median, ...
                                            obj.error.position_stats.std, ...
                                            obj.normstart_interval.time_stats.median,...
                                            1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(1),'bottom')

                ax(4) = subplot(2,3,5) ; hold on
                    h.errvel = plot(obj.normstart_interval.time,obj.error.velocity);
                    [hstd(2),~] = PlotPatch(obj.error.velocity_stats.median, ...
                                            obj.error.velocity_stats.std, ...
                                            obj.normstart_interval.time_stats.median,...
                                            1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(2),'bottom')

                ax(5) = subplot(2,3,3) ; hold on ; title('Integrated Error')
                    h.interrpos = plot(obj.normstart_interval.time,obj.int_error.position);
                    [hstd(3),~] = PlotPatch(obj.int_error.position_stats.median, ...
                                            obj.int_error.position_stats.std, ...
                                            obj.normstart_interval.time_stats.median,...
                                            1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(3),'bottom')

                ax(6) = subplot(2,3,6) ; hold on
                    h.interrvel = plot(obj.normstart_interval.time,obj.int_error.velocity);
                    [hstd(4),~] = PlotPatch(obj.int_error.velocity_stats.median, ...
                                            obj.int_error.velocity_stats.std, ...
                                            obj.normstart_interval.time_stats.median,...
                                            1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(4),'bottom')

                set(h.stimpos,      {'color'}, num2cell(obj.cmap,2))
                set(h.stimvel,      {'color'}, num2cell(obj.cmap,2))
                set(h.errpos,       {'color'}, num2cell(obj.cmap,2))
                set(h.errvel,       {'color'}, num2cell(obj.cmap,2))
                set(h.interrpos,    {'color'}, num2cell(obj.cmap,2))
                set(h.interrvel,    {'color'}, num2cell(obj.cmap,2))

                set([h.stimpos h.stimvel h.errpos h.errvel h.interrpos h.interrvel],'LineWidth',1)
                set(ax,'LineWidth',1,'FontWeight','bold')
                linkaxes(ax,'x')
                linkaxes(ax([1,3,6]),'y')
                %align_Ylabels(FIG)
            end
        end
        
    end
end

