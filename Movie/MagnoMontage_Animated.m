function [MOV] = MagnoMontage_Animated(rootpat,vidFs,body_ang,head_ang,rwing_ang,lwing_ang)
%% MagnoMontage_Animated: makes movie of animated fly in magnetic tether
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
rootpat = 'C:\Users\BC\Box\Git\Arena\Patterns';

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
FILE.montage    = [FILE.basename '_Animation.mp4'];
FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p'); % load DAQ data
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

FLY.nframe = length(TRIG.time_intrp_exp);
FLY.raw_center = [0 0];
radius = 150; % radius of pattern
thickness = 10; % radius display width

center = [0 0];
body_color = 'r';
head_color = 'c';
full_size = 150;
FlyAnimation = FlyModel(center, full_size, head_color, body_color, 1);

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    % VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'Uncompressed AVI');
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer',  'OpenGL', 'Units', 'inches','Position', [2 2 15 6]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
    
clear ax
ax(1) = subplot(6,8,[6:8,14:16,22:24,30:32,38:40,46:48]) ; cla; hold on % for fly animation
        axis image
ax(2) = subplot(6,8,[1:5,9:13]) ; cla ; hold on
        ylabel('Stimulus (°)','Color','w','FontSize',fontsize)
     	h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
ax(3) = subplot(6,8,[17:21,25:29]) ; cla ; hold on
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
ax(4) = subplot(6,8,[33:37,41:45]) ; cla ; hold on
        ylabel('Body (°)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.body = animatedline('Color','r','LineWidth',linewidth); % for body angle

set(ax(1),'Color',[204 255 255]./255)
set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(2:end-1), 'XTickLabel', [], 'XColor', 'none')
set(ax(2), 'YLim', 20*[-1 1], 'YTick', 15*[-1 0 1])
set(ax(3), 'YLim', 6 *[-1 1], 'YTick', 5 *[-1 0 1])
linkaxes(ax(2:end),'x')
align_Ylabels_ax(ax(2:end)')

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        disp(jj)
        
        % Display raw video
        subplot(6,8,[6:8,14:16,22:24,30:32,38:40,46:48]); cla; hold on ; axis off
            FlyAnimation = draw(FlyAnimation, FLY.int_body(jj), FLY.int_head(jj), ...
                                FLY.int_lwing(jj), FLY.int_rwing (jj), center);
                            
            % draw_semi_ellipse([0 0],50,0.75,0.9,FLY.int_body(jj),'m');
        
            % Make pattern ring
            make_pattern_ring(pattern_data.pattern,[PAT.pos_exp(jj),5],FLY.raw_center,...
                                    radius,thickness,'g');
    end
    
    % Pattern plot
 	subplot(6,8,[17:21,25:29]); hold on
        addpoints(h.pat, FLY.int_time(jj), PAT.norm(jj))
    
    % Head plot
    subplot(6,8,[17:21,25:29]); hold on
        addpoints(h.head, FLY.int_time(jj), FLY.int_head(jj))
    
   	% Body plot
    subplot(6,8,[33:37,41:45]); hold on
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