function [MOV] = montage_realtime_vid(rootdir,rootpat,vidFs,export)
%% montage_realtime: makes movie for fly in magnetic tether
%
% 	Includes fly video, registered video, body tracking, head tracking, 
%   wing tracking,& pattern position
%
%   INPUT:
%       rootdir     : directory containing all files
%       rootpat     : directory containing PATTERN file
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%
%   OUTPUT:
%       MOV         : structure containing movie 
%

clear ; clc ; close all
export = true;
% vidFs = 45;
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_body_reafferent\gain=0';
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
PATH.raw = rootdir;

if dirflag
    % Select data file
    [FILE.raw, PATH.head_track] = uigetfile({'*.mat'}, ...
        'Select fly file', PATH.raw, 'MultiSelect','off');
else
    FILE.raw = [mainfile , mainext];
end

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.montage = [FILE.basename '_montage.mp4'];

% Video data path
FILE.vid = fullfile(PATH.raw, [FILE.basename '.mp4']);

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Load data
disp('Loading Data ...')
pattern_data = load(fullfile(PATH.pat,FILE.pat),'pattern');
all_data = load(fullfile(PATH.raw,FILE.raw),'data','t_p','bAngles', 't_v', 'feedback');

% Load video data
Vreader = VideoReader(FILE.vid);
% data.vidData = read(Vreader, [1 Vreader.NumFrame]);

disp('DONE')

%% Get trigger times
daq_time = all_data.t_p;
daq_pattern = all_data.data(:,1);
trigger = all_data.data(:,3);
func_length = 60;
start_idx = nan;
n_frame = length(all_data.bAngles);
showplot = false;
[TRIG,PAT] = sync_pattern_trigger_realtime(daq_time, daq_pattern, func_length, ...
    trigger, start_idx, n_frame, true);

%% Get body angles offline
% [norm_ang,imgstats] = bodytracker(FLY.raw, 0, true, true, false);
% close all

%% Get video & body data
clc
FLY.frames = (TRIG.start_idx:TRIG.end_idx)';
FLY.n_frame = length(FLY.time);
FLY.time = TRIG.time;
[~,endI] = min(abs(FLY.time - func_length));
FLY.time_int = FLY.time(1:endI);
FLY.body = interp1(FLY.time, all_data.bAngles, FLY.time_int);
FLY.pat = interp1(PAT.time_sync, PAT.pos, FLY.time_int);
FLY.fs = 1 / mean(diff(FLY.time));
% FLY.body = all_data.bAngles(FLY.frames);
% FLY.pat = interp1(PAT.time_sync, PAT.pos, FLY.time);

% FLY.pat = PAT.pos_exp(FLY.frames);

% FLY.pat = interp1(PAT.time_sync, PAT.pos, FLY.time);

% FLY.body = FLY.body - median(FLY.body);
% FLY.pat = FLY.pat - median(FLY.pat);

close all
hold on
plot(FLY.time_int, FLY.pat, 'k', 'LineWidth', 1)
plot(FLY.time_int, FLY.body, 'r', 'LineWidth', 1)

%% For animation
pat_ypos = round(12*( median(all_data.data(:,3)) / 10));
off = 180;
pat_win = -40*[1 1];
radius = 108; % radius of pattern
thickness = 4; % radius display width
% FLY.body_centroid   = cat(1,imgstats(trig_range).Centroid); % centroid of image
% FLY.body_length     = cat(1,imgstats(trig_range).MajorAxisLength); % long axis of image
% FLY.body_tip        = FLY.body_centroid + repmat((FLY.body_length/4),1,2).* ...
%                                     [sind(FLY.int_body + off), -cosd(FLY.int_body + off)];
% FLY.reverse         = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
%                                      -[sind(FLY.int_body + off), -cosd(FLY.int_body + off)];
% FLY.heading = [FLY.body_centroid , FLY.body_tip];
% FLY.abdomen = [FLY.body_centroid , FLY.reverse];

%% Make Movie
clc
% Create video object
if export
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = FLY.fs;
    VID.Quality = 100;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 2000, 600]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
clear h ax ax_pat h_gaze
ax(1) = subplot(1,10,1:3) ; cla; hold on; axis image % for raw fly & pattern vid
        %title('Arena Frame','Color','w','FontSize', 1.2*fontsize)
        ax(1).XLim = [-pat_win(1) pat_win(1)+300];
        ax(1).YLim = [-pat_win(2) pat_win(2)+300];
        ax_pat(1) = axes; axis image
        set(ax_pat(1), 'Position', ax(1).Position, 'XLim', ax(1).XLim, 'YLim', ax(1).XLim)
        %ax(1).Position(3:4) = 0.9*ax(1).Position(3:4);
ax(2) = subplot(1,10,5:10)  ; cla ; hold on
        ylabel('(°)','Color','w','FontSize',fontsize)
        xlabel('Time (s)','Color','w','FontSize',fontsize)
     	h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
        h.body = animatedline('Color','r','LineWidth',linewidth); % for body angle
        leg = legend([h.body h.pat], 'body', 'display', 'Box', 'off', ...
            'TextColor', 'w', 'FontWeight', 'bold', 'FontSize', 14);
        ax(2).Position(4) = 0.7*ax(2).Position(4);
        ax(2).Position(2) = 2.2*ax(2).Position(2);
        %ax(2).XLim(1) = -5;
        
set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5, 'XLim', [-2 round(FLY.time_int(end) + 2)])
% set(ax(end),'XTick', 0:2:round(FLY.int_time(end)))

gs = pattern_data.pattern.gs_val + 1;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(ax_pat, cmap)

pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);
fig_sz = FIG.Position(3:4);

disp('Exporting Video...')
tic
frame = read(Vreader, 1);
set(FIG, 'CurrentAxes', ax(1)); cla; hold on; axis image
H.frame = imshow(frame);
H.body = [];
H.cent = [];
for n = 1:FLY.n_frame
    f = FLY.frames(n);
    if ~mod(n,100) || (n == 1)
        disp(n)
    end
    
    % Get frames
    rgb_frame = fliplr(read(Vreader, f));
    frame = rgb2gray(rgb_frame);
    frame = 1.1*imadjust(medfilt2(frame, [3 3])); % raw frame
    bw = imbinarize(frame);
    imgprop = regionprops(bw, 'Centroid', 'MajorAxisLength');
    FLY.body_centroid = cat(1,imgprop.Centroid); % centroid of image
    FLY.body_length = cat(1,imgprop.MajorAxisLength); % long axis of image
    FLY.body_tip = FLY.body_centroid + (FLY.body_length/4).* ...
                                        [sind(FLY.body(n) + off), -cosd(FLY.body(n) + off)];
    FLY.heading = [FLY.body_centroid , FLY.body_tip];
    pat_pos = FLY.pat(n);

    % Display raw video
    set(FIG, 'CurrentAxes', ax(1)); hold on; axis image
    delete(H.body)
    delete(H.cent)
        set(H.frame, 'CData', rgb_frame)
        H.body = plot([FLY.heading(1,1) FLY.heading(1,3)], ...
            [FLY.heading(1,2) FLY.heading(1,4)], 'r', 'LineWidth', 2);
        H.cent = plot(FLY.body_centroid(1), FLY.body_centroid(2), 'r.', 'MarkerSize', 25);

        if pat_ypos ~= 1
            % Make pattern ring
            set(FIG, 'CurrentAxes', ax_pat(1))
            cla(ax_pat(1))
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = radius * cos(theta) + 150;
            y = radius * sin(theta) + 150;
            z = zeros(1,length(x));
            hs = surface(ax_pat(1),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
        end

        if n == 1
            ax(1).Title.Position(2) = -20;
        end

    % Body & pattern plot
 	set(FIG, 'CurrentAxes', ax(2))
        addpoints(h.pat, FLY.time_int(n), FLY.pat(n))
        addpoints(h.body, FLY.time_int(n), FLY.body(n))
            
    if n == 1
        for y = 3:length(ax)
            ax(y).YLabel.Position(1) = -50;
        end
    end
    
    drawnow
    
    if export
        fig_frame = getframe(FIG);
        fig_frame.cdata = fig_frame.cdata(...
            round(0.1*fig_sz(2)):end-round(0.1*fig_sz(2)), ...
            round(0.12*fig_sz(1)):end-round(0.07*fig_sz(1)), :);
        writeVideo(VID,fig_frame);
    end
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')

%%
trig = all_data.data(:,3);
trig = 0.1*round(trig ./ 0.1);
trig(trig >= 3.1) = 3;
trig(trig < 2.85) = 0;
dxtrig = diff(trig);
dxtrig = [dxtrig(1) ; dxtrig];
[~,pks] = findpeaks(dxtrig, 'MinPeakHeight', 2);
pks_time = all_data.t_p(pks);
% pks_mag = trig(pks);
time_sync = all_data.t_p - pks_time(1);

% n_frame = length(bAngles);
% vid_time = linspace(pks_time(1), t_p(end), n_frame)' - pks_time(1);
vid_time = all_data.t_v;
pat = all_data.data(:,1);
pat = 3.75*round(96*(pat ./ 10));
pat = rad2deg(unwrap(deg2rad(pat)));

pat_pos = pat  - pat(1);
body = all_data.bAngles - mean(all_data.bAngles) + mean(pat_pos);

n_frame = 5500;
tintrp = vid_time(1:n_frame);
body_norm = interp1(vid_time, body, tintrp);
pat_norm = interp1(time_sync, pat_pos, tintrp);

ts = mean(diff(tintrp));
fs = 1 / ts;

fc = 10;
[b,a] = butter(5, fc/(fs/2), 'low');
body_filt = filtfilt(b, a, body_norm);
pat_filt = filtfilt(b, a, pat_norm);

pat_ypos = round(12*( median(all_data.data(:,2)) / 10));
end


