function [MOV] = MagnoMontage_SOS_v1(rootdir,rootpat,vidFs,export)
%% MagnoMontage_SOS_v1: makes movie for fly in magnetic tether
%
% 	Includes fly video, registered video, body tracking, head tracking, 
%   wing tracking,& pattern position
%
%   INPUT:
%       rootdir     : directory containing BENIFLY file
%       rootpat     : directory containing PATTERN file
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%
%   OUTPUT:
%       MOV         : structure containing movie 
%

% Example Input %
clear ; clc ; close all 
export = false;
vidFs = 50;
rootdir = 'E:\Experiment_SOS_v1';
rootpat = 'Q:\Box Sync\Git\Arena\Patterns';

% Create data paths
PATH.raw            = rootdir;
PATH.reg            = fullfile(PATH.raw,'registered');
PATH.body_track  	= fullfile(PATH.raw,'tracked_body');
PATH.head_track     = fullfile(PATH.reg,'tracked_head');
PATH.beninfly_track	= fullfile(PATH.reg,'tracked_head_wing');
PATH.mask           = fullfile(PATH.beninfly_track,'mask');

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Select tracked angle file (use head tracked files to select)
[FILE.raw, PATH.head_track] = uigetfile({'*.mat'}, ...
    'Select fly file', PATH.head_track, 'MultiSelect','off');

% Select pattern file
[FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
    'Select pattern file', rootpat, 'MultiSelect','off');

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.benifly   	= [FILE.basename '.csv'];
FILE.montage    = [FILE.basename '_Montage.mp4'];
FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p','vidData','t_v'); % load raw video & DAQ pattern positions
reg_data        = load(fullfile(PATH.reg,FILE.raw),'regvid','trf'); % load registered video
head_data    	= load(fullfile(PATH.head_track,FILE.raw),'hAngles','cPoint'); % load head angles
body_data    	= load(fullfile(PATH.body_track,FILE.raw),'bAngles','imgstats','initframe'); % load body angles
disp('DONE')

%% Get pattern data & sync with start of visual stimulus
pattern_total_time = 20;
[TRIG,PAT] = sync_pattern_trigger(raw_data.t_p, raw_data.data(:,2), pattern_total_time, ...
                        raw_data.data(:,1), true, false);
% Get kinematics data
FLY.time    = TRIG.time_sync; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 15; % cut off frequency for lpf
[b,a]       = butter(2,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.body    = body_data.bAngles; % body angles [deg]
FLY.head    = filtfilt(b,a,head_data.hAngles); % head angles [deg]
% FLY.head    = filtfilt(b,a,rad2deg(benifly_data.Head)); % head angles [deg]
FLY.head    = FLY.head - mean(FLY.head); % head angles [deg]
FLY.lwing   = rad2deg(hampel(FLY.time,benifly_data.LWing)); % left wing angles [deg]
FLY.rwing   = rad2deg(hampel(FLY.time,benifly_data.RWing)); % right wing angles [deg]
FLY.lwing   = filtfilt(b,a,FLY.lwing); % left wing angles [deg]
FLY.rwing   = filtfilt(b,a,FLY.rwing); % right wing angles [deg]
FLY.wba     = FLY.lwing - FLY.rwing; % delta wing-beat-amplitude [deg]
[b,a]       = butter(2,20/(FLY.Fs/2),'low'); % make lpf
FLY.wba     = filtfilt(b,a,FLY.wba); % delta wing-beat-amplitude [deg]
FLY.wba     = FLY.wba - mean(FLY.wba); % delta wing-beat-amplitude [deg]

% Normalize fly kinematics for experimental window
FLY.int_time    = TRIG.time_intrp_exp; % video time
FLY.int_body    = FLY.body(TRIG.range);
FLY.int_head    = FLY.head(TRIG.range);
FLY.int_lwing 	= FLY.lwing(TRIG.range);
FLY.int_rwing  	= FLY.rwing(TRIG.range);
FLY.int_wba     = FLY.wba(TRIG.range);
PAT.norm        = 3.75*(PAT.pos_exp - mean(PAT.pos_exp));

%% Get video data
FLY.raw = (flipud(squeeze(raw_data.vidData(:,:,TRIG.range)))); % raw video data
FLY.reg = squeeze(reg_data.regvid(:,:,TRIG.range)); % registered video data
FLY.raw_crop = FLY.raw;
FLY.reg_crop = FLY.reg;
FLY.nframe = size(FLY.raw,3);

[FLY.raw_xP,FLY.raw_yP,~] = size(FLY.raw_crop); % get size of raw video
[FLY.reg_xP,FLY.reg_yP,~] = size(FLY.reg_crop); % get size of registered video
FLY.raw_center = [round(FLY.raw_yP/2) , round(FLY.raw_xP/2)]; % center point for pattern & fly
FLY.reg_center = [round(FLY.reg_yP/2) , round(FLY.reg_xP/2)]; % center point for pattern & fly

radius = floor(max([FLY.raw_yP FLY.raw_xP])/1.8); % radius of pattern
thickness = 15; % radius display width

%% Get benifly parameters/mask
% Heading
fid = fopen(fullfile(PATH.mask, FILE.mask));
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
params = jsondecode(str);

FLY.wing_length = 150;

FLY.lwing_hinge = [params.gui.left.hinge.x  , params.gui.left.hinge.y];
FLY.rwing_hinge = [params.gui.right.hinge.x , params.gui.right.hinge.y];

FLY.lwing_tip = FLY.lwing_hinge - FLY.wing_length*[cosd(FLY.int_lwing),  sind(FLY.int_lwing)];
FLY.rwing_tip = FLY.rwing_hinge + FLY.wing_length*[cosd(FLY.int_rwing), -sind(FLY.int_rwing)];

FLY.head_length = 35;
FLY.head_hinge = [head_data.cPoint.X  , head_data.cPoint.Y];
% FLY.head_hinge = [params.gui.head.hinge.x  , params.gui.head.hinge.y];
FLY.head_tip   = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head) , -cosd(FLY.int_head)];

FLY.body_centroid   = cat(1,body_data.imgstats(TRIG.range).Centroid); % centroid of image
FLY.body_length     = cat(1,body_data.imgstats(TRIG.range).MajorAxisLength); % long axis of image
FLY.body_tip        = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
                                    [sind(FLY.int_body), -cosd(FLY.int_body)];
FLY.reverse         = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
                                     -[sind(FLY.int_body), -cosd(FLY.int_body)];
FLY.heading = [FLY.body_centroid , FLY.body_tip];
FLY.abdomen = [FLY.body_centroid , FLY.reverse];

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    % VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'Uncompressed AVI');
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    % VID.LosslessCompression = true;
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 1.75*16*40, 1.2*16*50]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
    
clear ax
ax(1) = subplot(14,2,1:2:15) ; cla; hold on; axis square % for raw fly & pattern vid
        title('Arena Frame','Color','w','FontSize',fontsize)
ax(2) = subplot(14,2,2:2:16) ; cla; hold on; axis square % for registered fly & pattern vid
        title('Body Frame','Color','w','FontSize',fontsize)
ax(3) = subplot(14,2,(16+1):(16+4))  ; cla ; hold on
        ylabel('Stimulus (°)','Color','w','FontSize',fontsize)
     	h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
ax(4) = subplot(14,2,(16+4+1):(16+2*4)) ; cla; hold on;
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
ax(5) = subplot(14,2,(16+2*4+1):(16+3*4)) ; cla; hold on;
        % ylabel('\Delta WBA (°)','Color','w','FontSize',fontsize)
        ylabel('Body (°)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.body = animatedline('Color','r','LineWidth',linewidth); % for body angle

set(ax(3:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(3:end-1), 'XTickLabel', [], 'XColor', 'none')
set(ax(3), 'YLim', 20*[-1 1], 'YTick', 15*[-1 0 1])
set(ax(4), 'YLim', 6 *[-1 1], 'YTick', 5 *[-1 0 1])
% set(ax(5), 'YLim', 40*[-1 1], 'YTick', 30*[-1 0 1])
linkaxes(ax(3:end),'x')
align_Ylabels_ax(ax(3:end)')

pp = 1;
iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        disp(jj)
        Frame.raw   = 2*FLY.raw(:,:,jj); % current raw frame        
        Frame.reg   = 2*FLY.reg(:,:,jj); % current registered frame
        
        % Display raw video
        subplot(14,2,1:2:15); cla; hold on; axis square
            imshow(Frame.raw)
           
            % Show ellipse
            % h.ellps = ellipse(FLY.body_centroid(jj,:), FLY.body_length(jj), 0.5, 0.90, FLY.int_body(jj), 'r');
            % delete([h.ellps{[1,3:6]}])
            % h.ellps{2}.Color = 'r';
            % h.ellps{2}.LineWidth = 1;

            h.heading = semi_ellipse(FLY.body_centroid(jj,:), FLY.body_length(jj)/2, 0.5, 0.90, ...
                                            180 - FLY.int_body(jj), 'r');
            delete([h.heading{2:4}])
            alpha(h.heading{1},0.3)
            h.heading{1}.LineStyle = 'none';
            
         	% Show centroid & heading
            % plot(FLY.heading(jj,[1,3]), FLY.heading(jj,[2,4]), '-r', 'LineWidth', 2) % centroid
            % plot(FLY.abdomen(jj,[1,3]), FLY.abdomen(jj,[2,4]), '-b', 'LineWidth', 2) % heading
        
            % Make pattern ring
            make_pattern_ring(pattern_data.pattern,[PAT.pos_exp(jj),5],FLY.raw_center,...
                                    radius,thickness,'g');
     	
        % Display registered video
        subplot(14,2,2:2:16); cla; hold on; axis square
            imshow(Frame.reg)            
          	plot([FLY.head_hinge(1), FLY.head_tip(jj,1)], ... % update line drawn to head
                 [FLY.head_hinge(2), FLY.head_tip(jj,2)], 'Color', 'c', 'LineWidth', 1.5)
                         
          	plot([FLY.lwing_hinge(1) , FLY.lwing_tip(jj,1)], ... % update line drawn to left wing
                [FLY.lwing_hinge(2) , FLY.lwing_tip(jj,2)], 'Color', 'm', 'LineWidth', 1.5)
            
        	plot([FLY.rwing_hinge(1) , FLY.rwing_tip(jj,1)], ... % update line drawn to right wing
                [FLY.rwing_hinge(2) , FLY.rwing_tip(jj,2)], 'Color', 'm', 'LineWidth', 1.5)
            
            plot(FLY.lwing_hinge(1),FLY.lwing_hinge(2),'m.','MarkerSize',10)
            plot(FLY.rwing_hinge(1),FLY.rwing_hinge(2),'m.','MarkerSize',10)
            
            % Make pattern ring
            make_pattern_ring(pattern_data.pattern,[PAT.pos_exp(jj),5],FLY.reg_center,...
                                    radius,thickness,'g');
    end
    
    % Pattern plot
 	subplot(14,2,(16+4+1):(16+2*4)); hold on
        addpoints(h.pat, FLY.int_time(jj), PAT.norm(jj))
    
    % Head plot
    subplot(14,2,(16+4+1):(16+2*4)); hold on
        addpoints(h.head, FLY.int_time(jj), FLY.int_head(jj))
    
   	% Body plot
    subplot(14,2,(16+2*4+1):(16+3*4)); hold on
    % addpoints(h.wba, FLY.int_time(jj), FLY.int_wba(jj))
        % addpoints(h.body, FLY.int_time(jj), FLY.int_body(jj))
        addpoints(h.body, FLY.int_time(jj), FLY.int_wba(jj))

    drawnow
    
    if expframe(jj)
        fig_frame = getframe(FIG);
        if export
            try
                writeVideo(VID,fig_frame);
            catch
                pause(1)
                writeVideo(VID,fig_frame);
            end
        end
        
     	% Store frame
        % MOV(pp) = fig_frame;
        
        pp = pp + 1;
    end
    pause(0.001)
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')
end