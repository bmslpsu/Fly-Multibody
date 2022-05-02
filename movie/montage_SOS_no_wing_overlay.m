function [MOV] = montage_SOS_no_wing_overlay(rootdir,rootpat,vidFs,export)
%% montage_SOS_no_wing_overlay: makes movie for fly in magnetic tether
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

% Example Input
clear ; clc ; close all
export = true;
vidFs = 50;
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2_head_fixed';
% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v2';
% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_vel_250';
% rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SS_amp_3.75';
rootpat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';

if ~isfolder(rootdir)
    dirflag = false;
    [rootdir,mainfile,mainext] = fileparts(rootdir);
else
    dirflag = true;
end

if ~isfolder(rootpat)
	[PATH.pat,FILE.pat,patext] = fileparts(rootpat);
    FILE.pat = [FILE.pat , patext];
else
	% Select pattern file
    [FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
        'Select pattern file', rootpat, 'MultiSelect','off');
end


% Create data paths
PATH.raw            = rootdir;
PATH.reg            = fullfile(PATH.raw,'registered');
PATH.func        	= fullfile(PATH.raw,'function');
PATH.body_track  	= fullfile(PATH.raw,'tracked_body');
PATH.head_track     = fullfile(PATH.reg,'tracked_head_tip');
PATH.beninfly_track	= fullfile(PATH.reg,'tracked_head_wing');
PATH.mask           = fullfile(PATH.beninfly_track,'');

if dirflag
    % Select tracked angle file (use head tracked files to select)
    [FILE.raw, PATH.head_track] = uigetfile({'*.mat'}, ...
        'Select fly file', PATH.head_track, 'MultiSelect','off');
    [FILE.func, PATH.func] = uigetfile({'*.mat'}, ...
        'Select function file', PATH.func, 'MultiSelect','off');
else
    FILE.raw = [mainfile , mainext];
    %FILE.func = [mainfile , mainext];
end

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.main       = [FILE.basename '.mat'];
FILE.benifly   	= [FILE.basename '.csv'];
FILE.montage    = [FILE.basename '_montage_no_wing_overlay.mp4'];
FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
func_data       = load(fullfile(PATH.func,FILE.func),'All'); % load function
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.main),'data','t_p','vidData','t_v'); % load raw video & DAQ pattern positions
reg_data        = load(fullfile(PATH.reg,FILE.main),'regvid','trf'); % load registered video
head_data    	= load(fullfile(PATH.head_track,FILE.main),'head_data','head_mask'); % load head angles
body_data    	= load(fullfile(PATH.body_track,FILE.main),'bAngles','imgstats','initframe'); % load body angles
disp('DONE')

%% Get pattern data & sync with start of visual stimulus
func_time = round(func_data.All.time(end));
startI = round(5000*0.5);
[TRIG,PAT] = sync_pattern_trigger(raw_data.t_p, raw_data.data(:,2), func_time, ...
                        raw_data.data(:,1), true, startI, false, false);
trig_time = TRIG.time_sync;

% Get kinematics data
FLY.time    = trig_time; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 20; % cut off frequency for lpf
[b,a]       = butter(3, FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.body    = body_data.bAngles; % body angles [°]
FLY.head    = head_data.head_data.angle; % head angles [°]
FLY.head    = filtfilt(b,a,FLY.head); % head angles [°]
FLY.lwing   = rad2deg(hampel(FLY.time,benifly_data.LWing)); % left wing angles [°]
FLY.rwing   = rad2deg(hampel(FLY.time,benifly_data.RWing)); % right wing angles [°]
FLY.lwing   = filtfilt(b,a,FLY.lwing); % left wing angles [°]
FLY.rwing   = filtfilt(b,a,FLY.rwing); % right wing angles [°]
FLY.wba     = FLY.lwing - FLY.rwing; % delta wing-beat-amplitude [°]
[b,a]       = butter(3, 20/(FLY.Fs/2), 'low'); % make lpf
FLY.wba     = filtfilt(b,a,FLY.wba); % delta wing-beat-amplitude [°]

% Normalize fly kinematics for experimental window
FLY.int_time    = TRIG.time_intrp_exp; % video time
FLY.int_body    = FLY.body(TRIG.range);
FLY.int_head    = FLY.head(TRIG.range);
FLY.int_gaze    = FLY.int_body + FLY.int_head;
FLY.int_lwing 	= FLY.lwing(TRIG.range);
FLY.int_rwing  	= FLY.rwing(TRIG.range);
FLY.int_wba     = FLY.wba(TRIG.range);

Fc_pat = 25;
[b,a] = butter(3, Fc_pat/(FLY.Fs/2), 'low'); % make lpf
pat_filt = filtfilt(b, a, 3.75*PAT.pos_intrp_exp);
PAT.norm = pat_filt - mean(pat_filt);

PAT.norm = interp1(func_data.All.time, func_data.All.X, FLY.int_time , 'pchip');

pat_ylim = 10*[floor(min(PAT.norm)./10) ceil(max(PAT.norm)./10)];
pat_ypos = round(12*( median(raw_data.data(:,3)) / 10));

%% Get video data
FLY.raw = fliplr(squeeze(raw_data.vidData(:,:,TRIG.range))); % raw video data
FLY.reg = squeeze(reg_data.regvid(:,:,TRIG.range)); % registered video data

[FLY.raw_xP,FLY.raw_yP,~] = size(FLY.raw); % get size of raw video
[FLY.reg_xP,FLY.reg_yP,~] = size(FLY.reg); % get size of registered video
FLY.raw_center = [round(FLY.raw_yP/2) , round(FLY.raw_xP/2)]; % center point for pattern & fly
FLY.nframe = size(FLY.raw,3);

% Make raw & reg same size
padsize = [FLY.raw_xP - FLY.reg_xP, FLY.raw_yP - FLY.reg_yP] ./ 2;
FLY.reg = padarray(FLY.reg, fix(padsize), 'pre');
FLY.reg = padarray(FLY.reg, ceil(padsize), 'post');

pat_win = -40*[1 1];
radius = floor(max([FLY.raw_yP FLY.raw_xP])/2.5); % radius of pattern
thickness = 6; % radius display width

%% Get benifly parameters/mask
fid = fopen(fullfile(PATH.mask, FILE.mask));
raw = fread(fid,inf); 
str = char(raw');
fclose(fid);
params = jsondecode(str);

FLY.body_reg = -90 - rad2deg(atan2(params.gui.head.hinge.y - params.gui.abdomen.hinge.y, ...
                            params.gui.head.hinge.x - params.gui.abdomen.hinge.x));

FLY.wing_length = 150;
FLY.lwing_hinge = fliplr(fix(padsize)) + [params.gui.left.hinge.x  , params.gui.left.hinge.y];
FLY.rwing_hinge = fliplr(fix(padsize)) + [params.gui.right.hinge.x , params.gui.right.hinge.y];
FLY.lwing_tip = FLY.lwing_hinge - FLY.wing_length*[cosd(FLY.int_lwing - FLY.body_reg),  ...
                                                    sind(FLY.int_lwing - FLY.body_reg)];
FLY.rwing_tip = FLY.rwing_hinge + FLY.wing_length*[cosd(FLY.int_rwing + FLY.body_reg), ...
                                                    -sind(FLY.int_rwing + FLY.body_reg)];

FLY.body_glob = head_data.head_mask.global;
FLY.head_length = 36;
FLY.head_hinge = fliplr(fix(padsize)) + head_data.head_mask.move_points.rot + [0 0];
FLY.head_tip   = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head + FLY.body_glob) , ...
                    -cosd(FLY.int_head + FLY.body_glob)];

FLY.body_centroid   = cat(1,body_data.imgstats(TRIG.range).Centroid); % centroid of image
FLY.body_length     = cat(1,body_data.imgstats(TRIG.range).MajorAxisLength); % long axis of image
FLY.body_tip        = FLY.body_centroid + repmat((FLY.body_length/4),1,2).* ...
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
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = vidFs;
    VID.Quality = 100;
    open(VID)
end

body_color = [1 0 0];
head_color = [0 0.4 1];
gaze_color = [0.5 0.3 1];

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.7*[100, 100, 1.75*16*40, 1.4*16*50]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
ph = 10;
rows = 9;
clear ax ax_pat
ax(1) = subplot(rows,2,1:2:(ph-1)) ; cla; hold on; axis image % for raw fly & pattern vid
        title('Arena Frame','Color','w','FontSize', 1.2*fontsize)
        ax(1).XLim = [-pat_win(1) pat_win(1)+FLY.raw_xP];
        ax(1).YLim = [-pat_win(2) pat_win(2)+FLY.raw_yP];
        H.raw = imshow(FLY.raw(:,:,1));
        H.heading = {};
      	H.body = [];
        colormap(ax(1), gray)
        freezeColors
        ax(1).Title.Position(2) = 0.9*ax(1).Title.Position(2);
        ax_pat(1) = axes; axis image
        set(ax_pat(1), 'Position', ax(1).Position, 'XLim', ax(1).XLim, 'YLim', ax(1).XLim, 'Color', 'none')
ax(2) = subplot(rows,2,2:2:ph) ; cla; hold on; axis image % for registered fly & pattern vid
        title('Body Frame','Color','w','FontSize', 1.2*fontsize)
        ax(2).XLim = [-pat_win(1) pat_win(1)+FLY.raw_xP];
        ax(2).YLim = [-pat_win(2) pat_win(2)+FLY.raw_yP];
        H.reg = imshow(FLY.reg(:,:,1));
        H.head = [];
        colormap(ax(2), gray)
        freezeColors
        ax(2).Title.Position(2) = ax(1).Title.Position(2);
        ax_pat(2) = axes; axis image
        set(ax_pat(2), 'Position', ax(2).Position, 'XLim', ax(2).XLim, 'YLim', ax(2).XLim, 'Color', 'none')
ax(3) = subplot(rows,2,(ph+1):(ph+4))  ; cla ; hold on
        ylabel('Stimulus & Body (°)','Color','w','FontSize',fontsize)
     	h.pat = animatedline('Color',[0.9 0.9 0.9],'LineWidth',linewidth); % for pattern angle
        h.body = animatedline('Color',body_color,'LineWidth',linewidth); % for body angle
        h.gaze = animatedline('Color', gaze_color,'LineWidth', linewidth); % for body angle
        h.body = animatedline('Color',body_color,'LineWidth',linewidth); % for body angle
ax(4) = subplot(rows,2,(ph+4+1):(ph+2*4)) ; cla; hold on
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color',head_color,'LineWidth',linewidth); % for head angle
        xlabel('Time (s)','Color','w','FontSize',fontsize)
        
set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
set(ax(3:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5, 'XLim', [-0.2 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.int_time(end)))
set(ax(3:end-1), 'XTickLabel', [], 'XColor', 'none')
% set(ax(5), 'YLim', 7*[-1 1], 'YTick', 5 *[-1 0 1])
head_ylim = 5*ceil(max(abs(FLY.int_head - median(FLY.int_head)))./5);
set(ax(4), 'YLim', (head_ylim)*[-1 1], 'YTick', head_ylim*[-1 0 1])
% set(ax(3), 'YLim', pat_ylim)
set(ax(3), 'YLim', 100*[-1 1])

linkaxes(ax(1:2), 'xy')
linkaxes(ax(3:end),'x')
align_Ylabels_ax(ax(3:end)')

% gs = pattern_data.pattern.gs_val + 1;
% cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
cmap = bone(2);
colormap(ax_pat(1), cmap)
colormap(ax_pat(2), cmap)

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);
bright_scale = 1.0;
fig_sz = FIG.Position(3:4);

disp('Exporting Video...')
tic
for jj = 1:FLY.nframe
    if expframe(jj)
        % Get frames
        %disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.raw = bright_scale*imadjust(median(FLY.raw(:,:,win), 3)); % raw frame
        Frame.reg = bright_scale*imadjust(median(FLY.reg(:,:,win), 3)); % registered frame
        pat_pos = 3.75*(mean(PAT.pos_intrp_exp(win)));
        
        % Display raw video
        set(FIG, 'CurrentAxes', ax(1)); hold on
            set(H.raw, 'CData', Frame.raw)
            % cellfun(@(x) delete(x), H.heading)
            % H.heading = semi_ellipse(mean(FLY.body_centroid(win,:),1), mean(FLY.body_length(win))/2, 0.5, 0.90, ...
            %                                 180 - mean(FLY.int_body(win)), 'r');
            % delete([H.heading{2:4}])
            % alpha(H.heading{1},0.3)
            % H.heading{1}.LineStyle = 'none';
            delete(H.body)
            H.body(1) = plot([mean(FLY.body_centroid(win,1),1) mean(FLY.body_tip(win,1),1)], ...
                 [mean(FLY.body_centroid(win,2),1) mean(FLY.body_tip(win,2),1)], ...
                        'Color', body_color, 'LineWidth', 2);
            H.body(2) = plot(mean(FLY.body_centroid(win,1),1), mean(FLY.body_centroid(win,2),1), '.', ...
                        'Color', body_color, 'MarkerSize', 15);
        
            if pat_ypos ~= 1
                % Make pattern ring
                set(FIG, 'CurrentAxes', ax_pat(1))
                cla(ax_pat(1))
                theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
                x = radius * cos(theta) + FLY.raw_center(1);
                y = radius * sin(theta) + FLY.raw_center(2);
                z = zeros(1,length(x));
                hs = surface(ax_pat(1),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                    'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
            end
     	   
        % Display registered video
        set(FIG, 'CurrentAxes', ax(2)); hold on
            set(H.reg, 'CData', Frame.reg)
            delete(H.head)
          	H.head = plot([FLY.head_hinge(1), mean(FLY.head_tip(win,1))], ... % update line drawn to head
                 [FLY.head_hinge(2), mean(FLY.head_tip(win,2))], 'Color', head_color, 'LineWidth', 2);
            plot(FLY.head_hinge(1), FLY.head_hinge(2), '.', 'Color', head_color, 'MarkerSize', 5)
            
            if pat_ypos ~= 1
                % Make pattern ring
                set(FIG, 'CurrentAxes', ax_pat(2)); cla
                hs = surface(ax_pat(2),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                    'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
            end
    end
    
    % Pattern & Body plot
 	set(FIG, 'CurrentAxes', ax(3))
        addpoints(h.pat, FLY.int_time(jj), PAT.norm(jj))
        %addpoints(h.gaze, FLY.int_time(jj), FLY.int_gaze(jj))
        addpoints(h.body, FLY.int_time(jj), FLY.int_body(jj) - median(FLY.int_body))

    % Head plot
    set(FIG, 'CurrentAxes', ax(4))
        addpoints(h.head, FLY.int_time(jj), 0*(FLY.int_head(jj) - median(FLY.int_head)))
        
    drawnow
    
    if jj == 1
        for y = 3:length(ax)
            ax(y).YLabel.Position(1) = -50;
        end
    end
    
    if export
        if expframe(jj)
            fig_frame = getframe(FIG);
            fig_frame.cdata = fig_frame.cdata(...
                round(0.07*fig_sz(2)):end-round(0.05*fig_sz(2)), ...
                round(0.03*fig_sz(1)):end-round(0.08*fig_sz(1)), :);
         	writeVideo(VID, fig_frame);
        end
    end
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')
end