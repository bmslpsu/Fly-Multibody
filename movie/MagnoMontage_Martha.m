function [MOV] = MagnoMontage_Martha(rootdir,vidFs,export)
%% MagnoMontage_SOS_v3: makes movie for fly in magnetic tether
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
vidFs = 160/3;
export = true;

rootdir = 'H:\EXPERIMENTS\MAGNO\TestMagno_Martha';

% Create data paths
PATH.raw            = rootdir;
PATH.reg            = fullfile(PATH.raw,'registered');
PATH.body_track  	= fullfile(PATH.raw,'tracked_body');
PATH.head_track     = fullfile(PATH.reg,'tracked_head');
%PATH.beninfly_track	= fullfile(PATH.reg,'tracked_head_wing');
%PATH.mask           = fullfile(PATH.beninfly_track,'mask');

[FILE.raw, PATH.head_track] = uigetfile({'*.mat'}, ...
    'Select fly file', PATH.head_track, 'MultiSelect','off');
    
% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
% FILE.benifly   	= [FILE.basename '.csv'];
FILE.montage    = [FILE.basename '_Montage.mp4'];
% FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
% benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p','vidData','t_v'); % load raw video & DAQ pattern positions
reg_data        = load(fullfile(PATH.reg,FILE.raw),'regvid','trf'); % load registered video
head_data    	= load(fullfile(PATH.head_track,FILE.raw),'hAngles','cPoint'); % load head angles
body_data    	= load(fullfile(PATH.body_track,FILE.raw),'bAngles','imgstats','initframe'); % load body angles
disp('DONE')

%% Get video data
FLY.raw = (flipud(squeeze(raw_data.vidData(:,:,:)))); % raw video data
FLY.reg = squeeze(reg_data.regvid(:,:,:)); % registered video data
FLY.raw_crop = FLY.raw;
FLY.reg_crop = FLY.reg;
FLY.nframe = size(FLY.raw,3);

[FLY.raw_xP,FLY.raw_yP,~] = size(FLY.raw_crop); % get size of raw video
[FLY.reg_xP,FLY.reg_yP,~] = size(FLY.reg_crop); % get size of registered video

%% Get pattern data & sync with start of visual stimulus
FLY.Fs = 160;
T = FLY.nframe * 1/FLY.Fs;

% Get kinematics data
FLY.time    = linspace(0,T,FLY.nframe); % video time
FLY.Fc      = 40; % cut off frequency for lpf
[b,a]       = butter(2,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.body    = body_data.bAngles; % body angles [deg]
FLY.head    = filtfilt(b,a,head_data.hAngles); % head angles [deg]
% FLY.head    = filtfilt(b,a,rad2deg(benifly_data.Head)); % head angles [deg]

% Normalize fly kinematics for experimental window
FLY.int_time    = FLY.time; % video time
FLY.int_body    = FLY.body;
FLY.int_head    = FLY.head;

%% Get benifly parameters/mask
% Heading
FLY.head_length = 65;
FLY.head_hinge = [head_data.cPoint.X  , head_data.cPoint.Y];
FLY.head_tip   = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head) , -cosd(FLY.int_head)];

FLY.body_centroid   = cat(1,body_data.imgstats(:).Centroid); % centroid of image
FLY.body_length     = cat(1,body_data.imgstats(:).MajorAxisLength); % long axis of image
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
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 1.75*16*40, 1.4*16*50]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
    
clear ax
ax(1) = subplot(4,2,[1,3]) ; cla; hold on; axis square % for raw fly & pattern vid
        title('Arena Frame','Color','w','FontSize',fontsize)
ax(2) = subplot(4,2,[2,4]) ; cla; hold on; axis square % for registered fly & pattern vid
        title('Body Frame','Color','w','FontSize',fontsize)
ax(3) = subplot(4,2,5:6) ; cla; hold on
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
ax(4) = subplot(4,2,7:8) ; cla; hold on
        % ylabel('\Delta WBA (°)','Color','w','FontSize',fontsize)
        ylabel('Body (°)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.body = animatedline('Color','r','LineWidth',linewidth); % for body angle

set(ax(3:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 (FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(3:end-1), 'XTickLabel', [], 'XColor', 'none')
set(ax(3), 'YLim', 11 *[-1 1], 'YTick', 5 *[-1 0 1])
% set(ax(4), 'YLim', 10 *[-1 1], 'YTick', 8 *[-1 0 1])
linkaxes(ax(3:end),'x')
% align_Ylabels_ax(ax(3:end)')

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.raw   = median(1*FLY.raw(:,:,win),3); % current raw frame        
        Frame.reg   = median(1*FLY.reg(:,:,win),3); % current registered frame
        
        % Display raw video
        subplot(4,2,[1,3]); cla; hold on; axis square
            imshow(Frame.raw)
            h.heading = semi_ellipse(median(FLY.body_centroid(win,:),1), median(FLY.body_length(win)/2), 0.5, 0.90, ...
                                            180 - median(FLY.int_body(win)), 'r');
            delete([h.heading{2:4}])
            alpha(h.heading{1},0.3)
            h.heading{1}.LineStyle = 'none';
     	
        % Display registered video
        subplot(4,2,[2,4]); cla; hold on; axis square
            imshow(Frame.reg)            
          	plot([FLY.head_hinge(1), FLY.head_tip(jj,1)], ... % update line drawn to head
                 [FLY.head_hinge(2), FLY.head_tip(jj,2)], 'Color', 'c', 'LineWidth', 1.5)
    end
      
    % Head plot
    subplot(4,2,5:6); hold on
        addpoints(h.head, FLY.int_time(jj), FLY.int_head(jj))
    
   	% Body plot
    subplot(4,2,7:8); hold on
        addpoints(h.body, FLY.int_time(jj), FLY.int_body(jj))
        
    drawnow
    
    if export
        if expframe(jj)
            fig_frame = getframe(FIG);
         	writeVideo(VID,fig_frame);
        end
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