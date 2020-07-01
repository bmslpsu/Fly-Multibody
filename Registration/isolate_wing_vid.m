function [wing_vid] = isolate_wing_vid(vid,showplot)
%% isolate_wing_vid: isolate and binarize wings
%
%   INPUT:
%       vid         : video data
%       showplot   	: show figures
%
%   OUTPUT:
%       wing_vid   	: isolated wing video data
%

vid = squeeze(vid);
dim = size(vid);
wing_vid = uint8(zeros(size(vid)));
SE = strel('disk',3,8);
edge_ratio = 0.2;
left_fill = edge_ratio*dim(2);
right_fill = dim(2) - left_fill;
for f = 1:dim(3)
    I = vid(:,:,f);
    IE = imdilate(I,SE);
    IE = medfilt2(IE,[10 10]);
    
    if f == 1
        h = imhist(IE);
        [t,~] = otsuthresh(h);
        thresh = 255*t;
        BW = IE;
        BW(BW >= thresh) = 255;
        BW(BW < thresh) = 0;
        stats = regionprops(BW,'BoundingBox','MinorAxisLength','MajorAxisLength');
        [~,bIidx] = max([stats.MajorAxisLength]);
        stats = stats(bIidx);
        
        up_cut = fix(1.5*stats.BoundingBox(2));
        low_cut = ceil(fix(stats.BoundingBox(2)) + stats.BoundingBox(4));
        left_cut = fix(stats.BoundingBox(1));
        right_cut = ceil(left_cut + stats.BoundingBox(3));
    end
    %left = IE(:,1:left_cut);
    %right = IE(:,right_cut:end);
    IE(:,left_cut:right_cut) = 0;
    BW2 = imbinarize(IE);
	BW2(:,1:left_fill) = true;
    BW2(:,right_fill:end) = true;
    BW2 = imfill(BW2, 'holes');
    BW2 = bwareaopen(BW2,50);
 	BW2 = 255*uint8(BW2);

    BW2(up_cut:low_cut,left_cut:right_cut) = I(up_cut:low_cut,left_cut:right_cut);
    
    if showplot
        figure (1)
        subplot(2,1,1) ; cla
            imshow(IE)
            rectangle('Position', stats.BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2)
        subplot(2,1,2) ; cla
            imshow(BW2)

        pause
    end
        
    wing_vid(:,:,f) = BW2;
end
    
end