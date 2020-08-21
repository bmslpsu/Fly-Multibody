
clear ; close all ; clc

root = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_forRoll';
% root = 'H:\EXPERIMENTS\RIGID\Experiment_SOS_v2';
[files, PATH] = uigetfile({'*.mat','files'}, 'Select files', root, 'MultiSelect','on');

load(fullfile(PATH, files),'vidData','t_v')

%%
obj = headtracker_area_v1(vidData);

%%
[~, ~, obj] = get_yaw(obj, 0.5);

%%
[obj] = iso_head(obj, 1);

%%
[~, obj] = get_roll(obj, 0.35, 36.33, true);

%% fly 10 trial 4
close all ; clc
[obj] = play_tracking(obj, 1, t_v);
% targetdir = fullfile(PATH, 'movie');
% [obj] = play_tracking(obj, 2, t_v, targetdir, files);

%%
clear A
test = uint8(zeros( size(obj.head_iso_vid,1), size(obj.head_iso_vid,2), 3, size(obj.head_iso_vid,3) ));
L = size(test,2);
cent = median( cat(1,obj.head_iso_props.Centroid), 1);
span = 0.15;
cut = floor(cent(1) - L/2*span) : ceil(cent(1) + L/2*span);
se_erode = strel('disk',1,8);
for n = 1:obj.dim(3)
    A = cell(1,1);
    A{1} = obj.head_iso_vid(:,:,n);
    
    mask = imbinarize(A{1});
    head = A{end}(mask);
    [thresh, em] = graythresh(head);
    
    A{end+1} = imbinarize(A{end}, 0.6235);
    A{end}(:,cut) = false;
    A{end+1} = bwareaopen(A{end}, 40);
    A{end+1} = imerode(A{end}, se_erode);
    A{end+1} = bwareaopen(A{end}, 40);
    
    stats = regionprops('table',A{end},'Centroid','PixelIdxList');
    
    [~,idx] = sort(stats.Centroid(:,1), 'ascend');
    
%     Is = length(idx);
%     if Is == 2
%         % pass
%         A{end+1} = A{end};
%     elseif Is < 2
%         warning('test')
%     elseif Is > 2
%     	frame = A{end};
%         for r = 2:Is-1
%             frame = frame;
%             mask = stats.PixelIdxList(idx(r));
%             frame(mask{1}) = uint8(0);
%         end
%        A{end+1} = frame;
%     end
%     A{end+1} = frame;
    %A{end+1}(A{end}) = 30;
    A{end+1} = imoverlay(A{1},A{end}, [0.6 0.1 0.1]);
    
    test(:,:,:,n) = A{end};
%     
%     montage(A)
%     pause(0.005)
end