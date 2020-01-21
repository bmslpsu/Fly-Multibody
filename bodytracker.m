function [] = bodytracker(vid)
%% bodytracker: calculates weights of redundant actuation in frequency domain
%
%   INPUT:
%       time   	:   time vector
%       ref     :	reference input
%       IOFv   	:	input-output frequncies present in data (leave empty to automatically detect)
%
%   OUTPUT:
%       MOV   	:   
%

vid = flipvid(vidData,'lr');
%%
[yp,xp,nframe] = size(vid);

body_angle = nan(nframe,1);
SE_erode = strel('disk',8,4);
SE_dilate = strel('disk',12,4);
offset = 0;
dthresh = 10;
tic
disp(1)
for idx = 1:nframe
    if ~mod(idx,100)
        disp(idx)
    end
    frame = vid(:,:,idx);
    bnframe = imbinarize(frame);
    bnframe = imerode(bnframe,SE_erode);
    bnframe = imdilate(bnframe,SE_dilate);
    [yfly,xfly] = find(bnframe==1);
    yfly = -yfly;
    
    figure (1); clf
    ax(1) = subplot(1,2,1); hold on ; cla ; axis square
        imshow(frame)
   	ax(2) = subplot(1,2,2); hold on ; cla ; axis square
        plot(xfly,yfly,'.')
        axis([1 xp -yp -1])

    ellipse_t = fit_ellipse(xfly,yfly,ax(2));
    body_angle(idx) = rad2deg(ellipse_t.phi) + offset;
    
%     if idx>1
%         if abs(body_angle(idx) - body_angle(idx-1)) > dthresh
%             if offset==0
%                 offset = 180;
%             elseif offset==180
%                 offset = 0;
%             else
%                 error('')
%             end
%         end
%     end
    
end
toc

%%
bnvid = logical(zeros(yp,xp,nframe));
tic
parfor idx = 1:nframe
    % idx
    frame = vid(:,:,idx);
    bnframe = imbinarize(frame);
    bnframe = imerode(bnframe,SE_erode);
    bnframe = imdilate(bnframe,SE_dilate);
    bnvid(:,:,idx) = logical(bnframe);
end
toc

%%
tic
for idx = 1:nframe
    if ~mod(idx,100)
        disp(idx)
    end
    bnframe = bnvid(:,:,idx);
    [yfly,xfly] = find(bnframe==1);
    yfly = -yfly;
    
    figure (1); clf
    ax(1) = subplot(1,2,1); hold on ; cla ; axis square
        imshow(bnframe)
   	ax(2) = subplot(1,2,2); hold on ; cla ; axis square
        plot(xfly,yfly,'.')
        axis([1 xp -yp -1])

    ellipse_t = fit_ellipse(xfly,yfly,ax(2));
    body_angle(idx) = rad2deg(ellipse_t.phi) + offset;
end
toc

end