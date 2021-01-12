function [TRIG,PAT] = sync_pattern_trigger(daq_time, daq_pattern, function_length, trigger, reg, start_idx, add1, debug)
%% sync_pattern_trigger: syncs camera frames with pattern for length of experiment
%
% 	Find where the pattern starts, find where frames are tiggered >>> align
% 	and sync times >>> bound range for length of experiment
%
%   INPUT:
%       daq_time        :   raw time recorded by DAQ
%       daq_pattern     :   pattern position in raw voltage (x or y channel) from DAQ
%       pattern_length  :   length of function [s]
%       trigger         :   camera trigger voltage signal from DAQ (rising edge is 
%                           trigger for each frame)
%       trig_center    	:   ratio of trigger width to use for frame times (0-1)
%       reg             :   BOOLEAN to set interpolated times for pattern (optional)
%       start_first     :   make pattern start at this index no matter what (optional)
%       add1          	:   BOOLEAN add 1st frame becuase missed the rising edge in the trigger signal
%       debug         	:   BOOLEAN show debug plot (optional)
%
%   OUTPUT:
%       TRIG          	: structure containing
%       PAT             : structure containing
%

if nargin < 8
    debug = false; % default
    if nargin < 7
        add1 = false; % default
        if nargin < 6
            start_idx = []; % default
            if nargin < 5
                reg = false; % default
            end
        end
    end
end

if isempty(reg)
   reg = false; % default
end

if isempty(add1)
   add1 = false; % default
end

% Camera trigger signal pulses and times
TRIG.pos                = round(trigger); % trigger values
TRIG.diff               = diff(TRIG.pos); % trigger derivative (rising edge triggers frame)
[TRIG.pks,TRIG.locs]    = findpeaks(TRIG.diff, 'MinPeakHeight', 1.5); % where each frame starts.
TRIG.locs               = TRIG.locs + 1; % this is where the frame starts
if add1
    TRIG.locs = [1 ; TRIG.locs]; % add 1st frame
    TRIG.pks  = [TRIG.pks(1) ; TRIG.pks]; % add 1st frame
end
TRIG.time = daq_time(TRIG.locs); % where each frame starts

% Convert pattern voltage to panel position, sync pattern start to 0 time,
% find where pattern ends
PAT.total_time          = function_length; % length of experiment
PAT.pos                 = round((96/10)*(daq_pattern)); % pattern position
PAT.diff                = diff(PAT.pos); % patten derivative

% Start of experiment
if isnan(start_idx)
    PAT.sync = TRIG.locs(1); % where first frame starts (start of experiment)
elseif ~isempty(start_idx)
    PAT.sync = start_idx; % set start point if specified
elseif isempty(start_idx)
    PAT.sync = find(abs(PAT.diff)>0,1,'first')+2; % where pattern first moves (start of experiment)
end

PAT.sync_time	= daq_time(PAT.sync); % start time
PAT.sync_val  	= PAT.pos(PAT.sync); % pattern value at sync
PAT.time_sync 	= daq_time - PAT.sync_time; % synced pattern time so 0 is the start
[~,PAT.end_idx]	= min(abs(PAT.time_sync - PAT.total_time)); % find where experiment ends

if PAT.end_idx >= length(daq_time)
    PAT.end_time = PAT.time_sync(PAT.end_idx+0); % end time
    %warning('Pattern to short')
else
    PAT.end_time = PAT.time_sync(PAT.end_idx+1); % end time
end

% Sync video with trigger & pattern
TRIG.time_sync          = TRIG.time - PAT.sync_time; % sync frames to start of visual motion
TRIG.sync_val           = TRIG.pos(PAT.sync); % trigger value at sync
TRIG.Fs                 = round( 1 / mean(diff(TRIG.time_sync )) ); % camera frame rate [Hz]

[~,TRIG.start_idx]    	= min(abs(TRIG.time_sync - 0)); % find where experiment starts
[~,TRIG.end_idx]    	= min(abs(TRIG.time_sync - PAT.total_time)); % find where experiment ends
TRIG.range              = TRIG.start_idx:TRIG.end_idx; % experiment range
TRIG.start_time        	= TRIG.time_sync(TRIG.start_idx); % start time
TRIG.end_time          	= TRIG.time_sync(TRIG.end_idx); % end time
TRIG.time_sync_exp     	= TRIG.time_sync(TRIG.start_idx:TRIG.end_idx); % frames during experiment.

PAT.pos_interp          = interp1(PAT.time_sync, PAT.pos, ... % interpolate pattern to match trigger
                                    TRIG.time_sync, 'nearest');
PAT.pos_exp             = PAT.pos_interp(TRIG.start_idx:TRIG.end_idx); % frames during experiment

if reg
    TRIG.time_intrp_exp = (0:(1/TRIG.Fs):PAT.total_time)'; % interpolated time
    PAT.pos_intrp_exp = interp1(PAT.time_sync, PAT.pos, ... % interpolate pattern to match regularized time
                               	TRIG.time_intrp_exp, 'nearest');       
end

% Debug sync
if debug
    fig = figure (100) ; clf
    set(fig,'Color','w','Units','inches','Position',[2 2 6 4])
    clear ax

    ax(1) = subplot(3,1,1) ; hold on ; axis tight ; title('Trigger')
        plot(PAT.time_sync,TRIG.pos,'k')
        plot(PAT.time_sync,[0;TRIG.diff],'b')
        if add1
            plot(TRIG.time_sync,TRIG.pks,'r.','MarkerSize',5)
        else
            plot(TRIG.time_sync,TRIG.pks,'r.','MarkerSize',5)
        end
        plot(0,TRIG.sync_val,'c.','MarkerSize',20)

        ylim([ax(1).YLim(1)-0.5 , ax(1).YLim(2)+0.5])

    ax(2) = subplot(3,1,2) ; hold on ; title('Sync Pattern Start to End')
        axis tight
        plot(daq_time,PAT.pos,'k','LineWidth',1)   
        plot(PAT.sync_time,PAT.sync_val,'c.','MarkerSize',20)

    ax(3) = subplot(3,1,3) ; hold on ; title('Sync Pattern End Time')
        axis tight
        plot(PAT.time_sync, PAT.pos, 'k', 'LineWidth', 1)
        plot(TRIG.time_sync, PAT.pos_interp, 'b', 'LineWidth', 1)
        plot(TRIG.time_sync_exp, PAT.pos_exp, 'g', 'LineWidth', 1)
        if reg
            plot(TRIG.time_intrp_exp, PAT.pos_intrp_exp, 'r', 'LineWidth', 1)
        end
        plot(PAT.end_time,PAT.pos(PAT.end_idx),'r.','MarkerSize',20)
        plot(0,PAT.sync_val,'c.','MarkerSize',20)

    linkaxes(ax(1:2),'x')
end
end
