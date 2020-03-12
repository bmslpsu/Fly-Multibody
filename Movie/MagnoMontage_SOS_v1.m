function [MOV] = MagnoMontage_SOS_v1(rootdir,rootpat,vidFs,export)
%% MagnoMontage_SOS_v1:  makes movie for fly in magnetic tether
% 	Includes fly video, head tracking, wing tracking, pattern position and plots of data
%   INPUT:
%       rootdir     : directory containing BENIFLY file
%       rootpat     : directory containing PATTERN file
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%   OUTPUT:
%       MOV         : structure containing movie 
%
% Example Input %
clear ; clc ; close all 
export = true;
vidFs = 50;
% rootdir = 'Q:\magno\Experiment_SOS';
% rootdir = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS';
rootdir = 'E:\Experiment_SOS_v1';
% rootpat = 'Q:\Box Sync\Git\Arena\Patterns';
rootpat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';

%%
% Create data paths
PATH.raw            = rootdir;
PATH.reg            = fullfile(PATH.raw,'registered');
PATH.body_track  	= fullfile(PATH.raw,'tracked_body');
PATH.head_track     = fullfile(PATH.reg,'tracked_head');
PATH.beninfly_track	= fullfile(PATH.reg,'tracked_headwing');

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Select tracked angle file (use head tracked files to select)
[FILE.raw, PATH.head_track] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select ANGLE file', PATH.head_track, 'MultiSelect','off');

% Select pattern file
[FILE.pat, PATH.pat] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select PATTERN file', rootpat, 'MultiSelect','off');

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.benifly   	= [FILE.basename '.csv'];
FILE.montage    = [FILE.basename '_Montage.avi'];

% Load data
fprintf('Loading Data... \n')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p','vidData','t_v'); % load raw video & DAQ pattern positions
reg_data        = load(fullfile(PATH.reg,FILE.raw),'regvid','trf'); % load registered video
head_data    	= load(fullfile(PATH.head_track,FILE.raw),'hAngles','cPoint'); % load head angles
body_data    	= load(fullfile(PATH.body_track,FILE.raw),'bAngles','imgstats'); % load body angles

fprintf('DONE \n')

%% Get pattern data % sync with start of visual stimulus
Pat.pos         = round((96/10)*(raw_data.data(:,2))); % pattern position
Pat.raw_time    = raw_data.t_p; % pattern time
Pat.diff        = diff(Pat.pos); % patten derivative
Pat.sync        = find(Pat.diff>0,1,'first')+1; % where pattern first moves (start of visual motion experiment)
Pat.sync_time   = Pat.raw_time(Pat.sync); % start time
Pat.sync_val    = Pat.pos(Pat.sync); % pattern value at sync
Pat.time        = Pat.raw_time - Pat.sync_time; % synced pattern time so 0 is the start
Pat.total_time  = 20; % length of experiment
[~,Pat.end_idx] = min(abs(Pat.time - Pat.total_time)); % find where experiment ends
Pat.end_time    = Pat.time(Pat.end_idx+1); % end time

% Sync video with trigger & pattern
Trig.raw_time   = raw_data.t_p; % DAQ raw times for trigger
Trig.pos        = round(raw_data.data(:,1)); % trigger values
Trig.diff       = diff(Trig.pos); % trigger derivative (rising edge triggers frame)
[pks,Trig.locs] = findpeaks(Trig.diff); % where each frame starts
Trig.locs       = [1 ; Trig.locs]; % add 1st frame
pks             = [pks(1) ; pks]; % add 1st frame
Trig.time       = Trig.raw_time(Trig.locs+1); % where each frame starts
Trig.time_sync  = Trig.time - Pat.sync_time; % sync frames to start of visual motion
% [~,Trig.start_idx] = min(abs(Trig.time_sync - 0)); % find where experiment ends
% [~,Trig.end_idx] = min(abs(Trig.time_sync - Pat.total_time)); % find where experiment ends
% Trig.time_sync 	= Trig.time_sync(Trig.start_idx:Trig.end_idx-1); % trigger times during experiment

% Get kinematics data
% Fly.time = Trig.time_sync; % video time
FLY.time    = Trig.time_sync; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 15; % cut off frequency for lpf
[b,a]       = butter(2,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.body    = body_data.bAngles;
FLY.head    = filtfilt(b,a,head_data.hAngles); % head angles [deg]
FLY.lwing   = rad2deg(hampel(FLY.time,benifly_data.LWing)); % left wing angles [deg]
FLY.rwing   = rad2deg(hampel(FLY.time,benifly_data.RWing)); % right wing angles [deg]
FLY.wba     = FLY.lwing - FLY.rwing; % delta wing-beat-amplitude [deg]
FLY.wba     = filtfilt(b,a,FLY.wba);
% FLY.wba     = FLY.wba - mean(FLY.wba);

% Normalize fly kinematics for experimental window
FLY.int_time    = (0:(1/FLY.Fs):Pat.total_time - (1/FLY.Fs))'; % video time
FLY.int_body    = interp1(FLY.time, FLY.body,  FLY.int_time, 'nearest'); % interpolate head to match fly video
FLY.int_head    = interp1(FLY.time, FLY.head,  FLY.int_time, 'nearest'); % interpolate head to match fly video
FLY.int_lwing   = interp1(FLY.time, FLY.lwing, FLY.int_time, 'nearest'); % interpolate left wing to match fly video
FLY.int_rwing   = interp1(FLY.time, FLY.rwing, FLY.int_time, 'nearest'); % interpolate head to match fly video
FLY.int_wba     = interp1(FLY.time, FLY.wba,   FLY.int_time, 'nearest'); % interpolate wba to match fly video

% Downsample pattern to match fly
Pat.int_pos = interp1(Pat.time, Pat.pos, FLY.int_time, 'nearest'); % interpolate pattern to match fly video
Pat.int_pos_deg = 3.75*(Pat.int_pos - mean(Pat.int_pos)); % interpolate pattern to match fly video

%% Get video data
FLY.vid_sync = find(FLY.time > 0,1,'first')-1;
FLY.vid_frames = round(FLY.vid_sync:(FLY.Fs*(Pat.total_time + (1./FLY.Fs)*FLY.vid_sync))-1)';
FLY.nframe = length(FLY.vid_frames);
FLY.raw = flipvid(squeeze(raw_data.vidData),'lr'); % raw video data
FLY.reg = squeeze(reg_data.regvid); % registered video data
FLY.raw_crop = FLY.raw;
FLY.reg_crop = FLY.reg;

% figure
% [FLY.raw_crop,FLY.raw_crop_area] = imcrop(FLY.raw(:,:,1,1));
% clf
% [FLY.reg_crop,FLY.reg_crop_area] = imcrop(FLY.reg(:,:,1,1));

% Debug sync
% fig = figure (100) ; clf
% set(fig,'Color','w','Units','inches','Position',[2 2 6 4])
% ax(1) = subplot(3,1,1) ; hold on
%     plot(Trig.raw_time,Trig.pos,'k')
%     plot(Trig.raw_time,[0;Trig.diff],'b')
%     plot(Trig.time,pks,'r*')
%     plot(Pat.raw_time,Pat.pos,'g')
%     plot(Pat.sync_time,Pat.sync_val,'c.','MarkerSize',20)
% ax(2) = subplot(3,1,2) ; hold on
%     plot(Pat.time,Pat.pos,'g')
% 	plot(Pat.end_time,Pat.pos(Pat.end_idx),'r*')
% ax(3) = subplot(3,1,3) ; hold on
%     plot(FLY.int_time,Pat.int_pos,'g')
%     plot(FLY.int_time,FLY.int_head,'b')
%     plot(FLY.int_time,FLY.int_lwing,'r')
% 	plot(FLY.int_time,FLY.int_rwing,'r')
% 	plot(FLY.int_time,FLY.int_wba,'k')

% [~,~,n_frame] = size(FLY.raw); % get size of video (assume raw & reg are same size)
[FLY.raw_xP,FLY.raw_yP,~] = size(FLY.raw_crop); % get size of raw video
[FLY.reg_xP,FLY.reg_yP,~] = size(FLY.reg_crop); % get size of registered video

FLY.raw_center = [round(FLY.raw_yP/2) , round(FLY.raw_xP/2)]; % center point for pattern & fly
FLY.reg_center = [round(FLY.reg_yP/2) , round(FLY.reg_xP/2)]; % center point for pattern & fly
FLY.head_R = 35;

radius = floor(max([FLY.raw_yP FLY.raw_xP])/1.8); % radius of pattern
thickness = 15; % radius display width.

%% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);
clc
% Create video object
if export
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'Uncompressed AVI');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 1.75*16*40, 1.2*16*50]);
linewidth = 1.25; % displayed line width
fontsize = 10;
    
clear ax
ax(1) = subplot(14,2,1:2:15) ; cla; hold on; axis square ; title('Arena Frame') % for raw fly & pattern vid
ax(2) = subplot(14,2,2:2:16) ; cla; hold on; axis square ; title('Body Frame') % for registered fly & pattern vid
ax(3) = subplot(14,2,(16+1):(16+4))  ; cla; hold on
        ylabel('Stimulus (°)','Color','w','FontSize',fontsize)
     	h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
ax(4) = subplot(14,2,(16+4+1):(16+2*4)) ; cla; hold on;
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
ax(5) = subplot(14,2,(16+2*4+1):(16+3*4)) ; cla; hold on;
        % ylabel('\Delta WBA (°)','Color','w','FontSize',fontsize)
        ylabel('Body (°)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.wba = animatedline('Color','r','LineWidth',linewidth); % for wba angle

set(ax(3:end), 'FontSize', 10, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(3:end-1), 'XTick', [])
set(ax(3), 'YLim', 20*[-1 1], 'YTick', 15*[-1 0 1])
set(ax(4), 'YLim', 6 *[-1 1], 'YTick', 5 *[-1 0 1])
% set(ax(5), 'YLim', 40*[-1 1], 'YTick', 30*[-1 0 1])
linkaxes(ax(3:end),'x')
align_Ylabels_ax(ax(3:end)')

pp = 1;
iter = round(FLY.Fs/vidFs); % how many frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==0,1); % logical to tell which frames to export
disp('Exporting Video...')
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        Frame.raw   = 2*FLY.raw(:,:,FLY.vid_frames(jj)); % current frame
        Disp.raw    = Frame.raw;
        % Disp.raw    = imcrop(Frame.raw,FLY.raw_crop_area); % raw video frame to display
        
        Frame.reg   = FLY.reg(:,:,FLY.vid_frames(jj)); % current frame
        Disp.reg    = 2*Frame.reg;
        % Disp.reg    = imcrop(Frame.reg,FLY.reg_crop_area); % registered video frame to display

        % Display raw video
        subplot(14,2,1:2:15); cla; hold on; axis square
            imshow(Disp.raw)
            
         	centroid = body_data.imgstats(jj).Centroid; % centroid of image
           	L = body_data.imgstats(jj).MajorAxisLength / 2; % long axis of image
          	head = centroid - L*[sind(FLY.int_body(jj)), cosd(FLY.int_body(jj))];
            heading = [centroid ; head];
            abdomen = centroid + L*[sind(FLY.int_body(jj)), cosd(FLY.int_body(jj))];
            reverse = [centroid ; abdomen];
            
            % Show ellipse
            h.ellps = ellipse(centroid, 2*L, 0.5, 0.90, FLY.int_body(jj), 'r');
            delete([h.ellps{[1,3:6]}])
            h.ellps{2}.Color = 'r';
            h.ellps{2}.LineWidth = 1;

            h.heading = semi_ellipse(centroid, -L, 0.5, 0.90, FLY.int_body(jj), 'r');
            delete([h.heading{2:4}])
            alpha(h.heading{1},0.3)
            h.heading{1}.LineStyle = 'none';
            
         	% Show centroid & heading
            plot(heading(:,1), heading(:,2), '-r', 'LineWidth', 2) % centroid
            plot(reverse(:,1), reverse(:,2), '-b', 'LineWidth', 2) % heading
            % plot(centroid(1),  centroid(2),  '.g', 'MarkerSize', 20) % reverse heading
        
        % Make pattern ring
        make_pattern_ring(pattern_data.pattern,[Pat.int_pos(jj),5],FLY.raw_center,...
                                    radius,thickness,'g');
     	
        % Display registered video
        subplot(14,2,2:2:16); cla; hold on; axis square
            imshow(Disp.reg)
            
            tracked.X = head_data.cPoint.X + FLY.head_R*sind(FLY.int_head(jj));
            tracked.Y = head_data.cPoint.Y  - FLY.head_R*cosd(FLY.int_head(jj));
            
          	plot([head_data.cPoint.X , tracked.X ], ... % update line drawn to head
                [head_data.cPoint.Y  , tracked.Y], 'Color', 'c', 'LineWidth', 1.5)
            
            
        % Make pattern ring
        make_pattern_ring(pattern_data.pattern,[Pat.int_pos(jj),5],FLY.reg_center,...
                                    radius,thickness,'g');

    end
    
    % Pattern plot
 	subplot(14,2,(16+4+1):(16+2*4)); hold on
        addpoints(h.pat, FLY.int_time(jj), Pat.int_pos_deg(jj))
    
    % Head plot
    subplot(14,2,(16+4+1):(16+2*4)); hold on
        addpoints(h.head, FLY.int_time(jj), FLY.int_head(jj))
    
   	% Wing plot
    subplot(14,2,(16+2*4+1):(16+3*4)); hold on
    % addpoints(h.wba, FLY.int_time(jj), FLY.int_wba(jj))
        addpoints(h.wba, FLY.int_time(jj), FLY.int_body(jj))

    drawnow
    
    if expframe(jj)
        if export
            writeVideo(VID,getframe(FIG));
        end
        
     	% Store frame
        MOV(pp) = getframe(FIG);
        pp = pp + 1;
    end
    pause(0.01)
end

if export
%  	disp('Saving...')
    close(VID) % close .avi
    Fs = FLY.Fs;
% 	save([root.mov dirName '.mat'],'MOV','Fs','-v7.3','-nocompression') % save movie as .mat file
end
disp('DONE')
end