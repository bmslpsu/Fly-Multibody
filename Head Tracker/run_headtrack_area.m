
clear ; close all ; clc

root = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_forRoll';
[files, PATH] = uigetfile({'*.mat','files'}, 'Select files', root, 'MultiSelect','on');

load(fullfile(PATH, files),'vidData','t_v')

%%
obj = headtracker_area_v1(vidData);

%%
[~, ~, obj] = get_yaw(obj, 0.5);

%%
[~, obj] = get_roll(obj, 0.35, 36.33);

%% fly 10 trial 4
close all ; clc
% [obj] = play_tracking(obj, 1, t_v);
targetdir = fullfile(PATH, 'movie');
[obj] = play_tracking(obj, 2, t_v, targetdir, files);

%%
clear A
test = obj.head_iso_vid(:,:,n);
for n = 1:obj.dim(3)
    A = cell(1,1);
    A{1} = 1.2*obj.head_iso_vid(:,:,n);
    
    mask = imbinarize(A{1});
    head = A{end}(mask);
    [thresh, em] = graythresh(head);
    
    A{end+1} = imbinarize(A{end}, 0.8);
    A{end+1} = bwareaopen(A{end}, 40);
    
    stats = regionprops('table',A{end},'Centroid','PixelIdxList');
    
    [~,idx] = sort(stats.Centroid(:,1), 'ascend');
    
    Is = length(idx);
    if Is == 2
        % pass
        A{end+1} = A{end};
    elseif Is < 2
        warning('test')
    elseif Is > 2
    	frame = A{end};
        for r = 2:Is-1
            frame = frame;
            mask = stats.PixelIdxList(idx(r));
            frame(mask{1}) = uint8(0);
        end
       A{end+1} = frame;
    end
    %A{end+1} = frame;
    A{1}(A{end}) = 30;
    
    test(:,:,n) = A{1};
    
%     montage(A)
%     pause(0.005)
end