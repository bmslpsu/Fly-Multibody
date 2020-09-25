function [heading,flip,out_frame,bias,ax,fig] = findheading(frame, debug, ratio)
%% findheading: find heading of fly in grey-scale image
%
%   INPUT:
%       frame    	:   frame to extract heading
%       debug       :   show debug figure 0=never, 1=always, 2=if close call (default = 2)
%      	ratio       :   frame ratio to compare top & bottom (default = 1/4)
%
%   OUTPUT:
%       heading  	:   heading orientation angle [°]
%       flip      	:   should we flip the heading? (0 = no: head at top, 
%                       1 = yes: head at bottom)
%       bias        :   ratio between upper and lower quadrant (values
%                       close to 1 are more ambiguous)
%       ax          :   axes handles
%       fig         :   figure handle
%

if nargin < 3
    ratio = 1/4; % default
    if nargin < 2
        debug = 2; % default
    end
end

SE_erode = strel('disk',8,8); % erosion mask

bnframe = imbinarize(frame); % binarize
bnframe = imerode(bnframe,  SE_erode); % erode

% Get image reigon stats
imgstats = regionprops(bnframe,'BoundingBox','Orientation','Image'); % image reigon properties
heading = imgstats(1).Orientation;

% Extract bounding reigon and rotate to 90°
check_frame = imrotate(imgstats(1).Image, 90 - imgstats(1).Orientation, 'loose');
head_frame = imrotate(frame, 90 - imgstats(1).Orientation, 'crop');
out_frame = head_frame;

% Get image reigon stats and extract bounding reigon of rotated reigon
imgstats = regionprops(check_frame,'Centroid','Orientation','Image'); % image reigon properties
[~,main_image_idx] = max(cellfun(@numel,{imgstats.Image}));
check_frame = imgstats(main_image_idx).Image;

imgstats = regionprops(imbinarize(head_frame),'BoundingBox','Image'); % image reigon properties
[~,flyIdx] = max(cellfun(@(x) numel(x), {imgstats.Image}));
head_frame = imcrop(head_frame,imgstats(flyIdx).BoundingBox);

% Get size & centroid of image
dim = size(check_frame);
imgstats = regionprops(check_frame,'Centroid'); % image reigon properties

% Make masks for the top & bottom 1/4 of image
mask_top = false(size(check_frame));
mask_bot = false(size(check_frame));
mid = fix( dim(1) / 2 );
top = fix( ratio * dim(1) );
bot = fix( (1-ratio) * dim(1) );
mask_top(1:top,:) = true;
mask_bot(bot:dim(1),:) = true;
% mask_full = mask_top | mask_bot;

% Get image reigons using masks
top_frame  = check_frame;
bot_frame  = check_frame;
% full_frame = check_frame;

top_frame(~mask_top)  = false;
bot_frame(~mask_bot)  = false;
% full_frame(~mask_full) = false;

% Calculate the bias between upper & lower reigons >>> compute the
% heading orientation 
bias = sum(top_frame) / sum(bot_frame);
if bias < 1 % head is at the top (correct guess)
    flip = false;
    top_color = 'r';
    bot_color = 'b';
    new_frame = head_frame;
elseif bias > 1 % head is at the bottom (inccorrect guess)
    flip = true;
    % Compute ajusted heading
    heading = heading + 180;
    if heading < 0 % keep angle positive
        heading = heading + 360;
    end
    top_color = 'b';
    bot_color = 'r';
    
    new_frame = flipud(head_frame); % flip frame so head is on top
    out_frame = flipud(out_frame); % flip frame so head is on top
end

% Check for ambiguous estimates % bring up debug window if heading estimate
% is too close to call
if debug == 2
    if abs(1-bias) < 0.1
        warning('Initial heading estimate may be incorect')
        debug = true;
    else
        debug = false;
    end
end

% Show plots
if debug
    fig(1) = figure (122); clf
    set(fig, 'Color', 'w', 'Units', 'inches', 'Name', 'Heading Detection')
    fig.Position(1:4) = [2 2 4 5];
        ax(1) = subplot(1,2,1); cla; hold on; axis image ; title(['Bias = ' num2str(bias)])
            imshow(check_frame) 
            plot(imgstats.Centroid(1),imgstats.Centroid(2),'.c','MarkerSize',30)
            plot([0 size(check_frame,2)], [imgstats.Centroid(2) imgstats.Centroid(2)], 'c', 'LineWidth', 2)
            plot([0 size(check_frame,2)], [mid mid], 'g', 'LineWidth', 2)
            plot([0 size(check_frame,2)], [top top], 'Color', top_color, 'LineWidth', 2)
            plot([0 size(check_frame,2)], [bot bot], 'Color', bot_color, 'LineWidth', 2)
            hold off

        ax(2) = subplot(1,2,2); cla ; hold on; axis image ; title('Adjusted Heading')
            imshow(new_frame)
            hold off
        
	set(ax,'FontSize',10)
    
	uicontrol('Style','pushbutton','String','Flip',    'Callback',@flipcheck,'Position', [75 50 60 20])
 	uicontrol('Style','pushbutton','String','Continue','Callback',@skipcheck,'Position', [250 50 60 20])
   	waitfor(fig)
else
    fig = [];
    ax  = [];
end

% Callback functions for flip & continue buttons
%------------------------------------------------
function flipcheck(~,~)
    close(fig)
    flip = ~flip; % invert the flip estimate
end

function skipcheck(~,~)
    close(fig) % continue
end
end
