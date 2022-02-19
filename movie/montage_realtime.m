function [MOV] = montage_realtime(rootdir,rootpat,vidFs,export)
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
vidFs = 50;
pat_ypos = 5;
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_reafferent';
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

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.montage    = [FILE.basename '_montage.mp4'];

% Load data
disp('Loading Data ...')
pattern_data = load(fullfile(PATH.pat,FILE.pat),'pattern');
data = load(fullfile(PATH.raw,FILE.raw),'data','t_p','vidData','t_v','bAngles', 'feedback');
disp('DONE')

%% Sync video & display data & filter
fly_time = data.t_v;
daq_time = data.t_p;

body = data.bAngles;
display = data.data(:,1);
display = 3.75*(round(96*(display ./ 10)) - 1);
display = rad2deg(unwrap(deg2rad(display)));

dx_pat = diff(display);
syncI = find(abs(dx_pat) > 5, 1, 'first') + 1;
sync_time = daq_time(syncI);

daq_time_sync = daq_time - sync_time;
display_intrp = interp1(daq_time_sync, display, fly_time, 'nearest');
% body_intrp = interp1(daq_time_sync, body, fly_time, 'nearest');

% Filter
Fs = 1 / mean(diff(fly_time));
Fc_low = 3;
[b_low,a_low] = butter(5, Fc_low/(Fs/2), 'low');
% [b_high,a_high] = butter(5, 0.3/(Fs/2), 'high');
display_filt = filtfilt(b_low, a_low, display_intrp);
body_filt = filtfilt(b, a, body);

%% Make perturbation function
func = 37.5*sin(2*pi*0.5*fly_time);




%% Get kinematics data
FLY.time = trig_time; % video time
FLY.Fs = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.body = data.bAngles - data.bAngles(1); % body angles [°]

% Fc_pat = 12;
% [b,a] = butter(3, Fc_pat/(FLY.Fs/2), 'low'); % make lpf
% PAT.norm = filtfilt(b, a, PAT.norm);

% Normalize fly kinematics for experimental window
trig_range = 1:length(FLY.body);
FLY.int_time = FLY.time(trig_range);
FLY.int_time = FLY.int_time - FLY.int_time(1);
FLY.int_body = FLY.body(trig_range);

PAT.int_norm = PAT.norm(trig_range);
PAT.int_norm = PAT.int_norm - PAT.int_norm(1);

%% Get video data
FLY.raw = fliplr(squeeze(data.vidData(:,:,trig_range))); % raw video data
[FLY.raw_xP,FLY.raw_yP,~] = size(FLY.raw); % get size of raw video
FLY.raw_center = [round(FLY.raw_yP/2) , round(FLY.raw_xP/2)]; % center point for pattern & fly
FLY.nframe = size(FLY.raw,3);

pat_win = -40*[1 1];
radius = floor(max([FLY.raw_yP FLY.raw_xP])/2.8); % radius of pattern
thickness = 4; % radius display width

%% Get body angles offline
[norm_ang,imgstats] = bodytracker(FLY.raw, 0, true, true, false);
close all

%% For animation
off = 0;
FLY.body_centroid   = cat(1,imgstats(trig_range).Centroid); % centroid of image
FLY.body_length     = cat(1,imgstats(trig_range).MajorAxisLength); % long axis of image
FLY.body_tip        = FLY.body_centroid + repmat((FLY.body_length/4),1,2).* ...
                                    [sind(FLY.int_body + off), -cosd(FLY.int_body + off)];
FLY.reverse         = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
                                     -[sind(FLY.int_body + off), -cosd(FLY.int_body + off)];
FLY.heading = [FLY.body_centroid , FLY.body_tip];
FLY.abdomen = [FLY.body_centroid , FLY.reverse];

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 2000, 600]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;
clear h ax ax_pat h_gaze
ax(1) = subplot(1,10,1:3) ; cla; hold on; axis square % for raw fly & pattern vid
        %title('Arena Frame','Color','w','FontSize', 1.2*fontsize)
        ax(1).XLim = [-pat_win(1) pat_win(1)+FLY.raw_xP];
        ax(1).YLim = [-pat_win(2) pat_win(2)+FLY.raw_yP];
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
    'LineWidth', 1.5, 'XLim', [-2 round(FLY.int_time(end))])
% set(ax(end),'XTick', 0:2:round(FLY.int_time(end)))

gs = pattern_data.pattern.gs_val + 1;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(cmap)

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);

fig_sz = FIG.Position(3:4);

disp('Exporting Video...')
tic
for jj = 1:FLY.nframe
    if ~mod(jj,100) || (jj == 1)
        disp(jj)
    end
    if expframe(jj)
        % Get frames
        %disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.raw = 1.1*imadjust(medfilt2(median(FLY.raw(:,:,win), 3), [3 3])); % raw frame
        pat_pos = (mean(PAT.int_norm(win)));
        
        % Display raw video
        set(FIG, 'CurrentAxes', ax(1)); cla; hold on; axis image
            imshow(Frame.raw)
            plot([mean(FLY.heading(win,1)) mean(FLY.heading(win,3))], ...
                [mean(FLY.heading(win,2)) mean(FLY.heading(win,4))], 'r', 'LineWidth', 2)
            plot(mean(FLY.body_centroid(win,1)), mean(FLY.body_centroid(win,2)), 'r.', 'MarkerSize', 25)
            %plot(mean(FLY.body_tip(win,1)), mean(FLY.body_tip(win,2)), 'c.', 'MarkerSize', 15)
        
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
     	
            if jj == find(expframe,1,'first')
                ax(1).Title.Position(2) = -20;
            end
    end
    
    % Body & pattern plot
 	set(FIG, 'CurrentAxes', ax(2))
        addpoints(h.pat, FLY.int_time(jj), PAT.int_norm(jj))
        addpoints(h.body, FLY.int_time(jj), FLY.int_body(jj))
        
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
                round(0.1*fig_sz(2)):end-round(0.1*fig_sz(2)), ...
                round(0.12*fig_sz(1)):end-round(0.07*fig_sz(1)), :);
         	writeVideo(VID,fig_frame);
        end
    end
    %pause(0.001)
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')
end