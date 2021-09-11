function [] = process_bizzi_head_eye()
%% process_bizzi_head_eye:
%
clear
close all ; clc

root = 'C:\Users\boc5244\OneDrive - The Pennsylvania State University\Research\Manuscripts\Multibody\figure\v4\';
image_file = fullfile(root, 'Bizzi_1973.PNG');
config_file = 'bizzi_config.mat';
config_path = fullfile(root, config_file);
config_check = exist(config_path, 'file');
data_path = fullfile(root, 'bizzi_data.mat');

% Read image
G = imread(image_file);
G = rgb2gray(G);
B{1} = G;
% B{end+1} = medfilt2(B{end}, [3 3]);
B{end+1} = imbinarize(B{end});
% montage(B)
P = B{end};

% ROI reigons
remove_rect_names = ["labels_right", "labels_left", "x_scale", "y_scale", "overlap_data"];
n_rect = length(remove_rect_names);
cc = hsv(n_rect);

global H
if config_check ~= 2 % set congig manually
    fig = figure (109);
    imshow(P)
    uicontrol('Style', 'pushbutton', 'String', 'Done', 'CallBack', @buttonPushed);

    % Draw ROI's
    rect = [];
    for n = 1:n_rect
        title(remove_rect_names(n), 'Interpreter', 'none')
        rect(n) = drawrectangle(gca, 'Color', cc(n,:), 'FaceAlpha', 0.15, 'LineWidth', 1, ...
            'Label', remove_rect_names(n));
        if n == n_rect
           title('Click done ROI''s are set')
        end
    end
    uiwait(fig) % wait to set ROI's
    save(config_path, 'H') % save when ROI's are set
else % load config
    load(config_path, 'H')
end

% Show ROI's
imshow(P)
% rect = [];
for n = 1:n_rect
    rect(n) = drawrectangle(gca, 'Color', cc(n,:), 'FaceAlpha', 0.15, 'LineWidth', 1, ...
        'Label', remove_rect_names(n), 'Position', H.(remove_rect_names(n)));
end

% Get scale images and calibrate
I.x_scale = imcrop(P, H.x_scale);
I.y_scale = imcrop(P, H.y_scale);
I.x_scale_sum = imcomplement(all(I.x_scale,1));
I.y_scale_sum = imcomplement(all(I.y_scale,2));

% Get calibration normalization based on axis scales
Data = [];
Data.x_cal = 0.1;
Data.y_cal = 40;
Data.x_scale = find(I.x_scale_sum, 1, 'last') - find(I.x_scale_sum, 1, 'first');
Data.y_scale = find(I.y_scale_sum, 1, 'last') - find(I.y_scale_sum, 1, 'first');
Data.x_norm = Data.x_cal ./ Data.x_scale;
Data.y_norm = Data.y_cal ./ Data.y_scale;

% Get data images
I.data = imcomplement(P);
for n = 1:n_rect-1
    mask = createMask(rect(n));
    I.data(mask) = false;
    %I.data = imerase(H.(remove_rect_names(n)));
end
sep_head_free = round(H.overlap_data(1));
sep_head_fixed = round([H.overlap_data(1) + H.overlap_data(3), ...
                        H.overlap_data(2) + H.overlap_data(4), ...
                        H.overlap_data(1)]);
I.data = bwareaopen(I.data, 5);

% Get head-free data image
I.head_free = I.data;
mask = createMask(rect(n_rect));
I.head_free(mask) = false;
I.head_free(:,1:sep_head_free) = false;

I.head_fixed = I.data;
I.head_fixed(:,sep_head_fixed(1):end) = false;
I.head_fixed(sep_head_fixed(2):end,sep_head_fixed(3):end) = false;

I.head_free = imclose(I.head_free, 10);
I.head_fixed = imclose(I.head_fixed, 10);

% Get raw data
[n_y, n_x] = size(I.data);
Data.head_free.x = (1:n_x)'; % x-data

[y_med_sort] = cluster_data(I.head_free, 3); % head free [head , eyes , gaze]
head2gaze = [1125 , 282]; % eyes transition to gaze in clusters at this point
y_med_class = y_med_sort;
y_med_class(head2gaze(1):end,2) = y_med_sort(head2gaze(1):end,3);
y_med_class(head2gaze(1):end,3) = y_med_sort(head2gaze(1):end,2);
Data.head_free.y = y_med_class;

[y_med_sort] = cluster_data(I.head_fixed, 1); % head fixed
Data.head_fixed.y = y_med_sort;

%% Normalize data
nanI = all(isnan(Data.head_free.y),2);
head_free_norm = Data.head_free.y(~nanI,:);
head_free_time = Data.head_free.x(~nanI);
head_free_time = head_free_time - head_free_time(1);
norm_yx = n_y - median(head_free_norm(1:10,:), 'all');

Data.head_free.time = (head_free_time - head_free_time(1)) * Data.x_norm;
Data.head_free.ts = mean(diff(Data.head_free.time));
Data.head_free.head = (n_y - head_free_norm(:,3) - norm_yx) * Data.y_norm;
Data.head_free.eyes = (n_y - head_free_norm(:,2) - norm_yx) * Data.y_norm;
Data.head_free.gaze = (n_y - head_free_norm(:,1) - norm_yx) * Data.y_norm;
Data.head_free.input = Data.head_free.gaze(end)*ones(size(Data.head_free.time));
Data.head_free.error = Data.head_free.input - Data.head_free.gaze;
%Data.head_fixed.eyes = (n_y - head_free_norm(:,1) - norm_yx) * Data.y_norm;
save(data_path, 'Data')

%% Plot
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 3])
clear h cc
cc.head = [0.9 0 0];
cc.eyes = [0 0.4 1];
cc.gaze = [0.5 0.3 1];
ax = subplot(1,1,1); cla ; hold on
    h(1) = plot(Data.head_free.time, Data.head_free.head, 'Color', cc.head);
    h(2) = plot(Data.head_free.time, Data.head_free.eyes, 'Color', cc.eyes);
    h(3) = plot(Data.head_free.time, Data.head_free.gaze, 'Color', cc.gaze);
    h(4) = plot(Data.head_free.time, Data.head_free.input, 'Color', 'k');
    h(5) = plot(Data.head_free.time, Data.head_free.error, 'Color', 'g');

set(h, 'LineWidth', 1)
set(ax, 'Color', 'w', 'LineWidth', 0.75, 'Box', 'on')
set(ax, 'XGrid', 'off', 'YGrid', 'off')
set(ax, 'YLim', [-2 30])
set(ax, 'XLim', [-0.01 0.4])

xlabel('time (s)')
ylabel('(°)')

%% System ID
u = Data.head_free.input;
e = Data.head_free.error;
y = Data.head_free.head;
ts = Data.head_free.ts;
data = iddata(y, u, ts);
% d = delayest(data, 1, 1, 0, 0.1)
sys = tfest(data, 2, 2, nan);

step(sys)

%% Debug
close all
subplot(1,1,1) ; cla
    imshow(I.data) ; hold on
    plot(Data.head_free.x, Data.head_free.y, '-', 'MarkerSize', 10, 'LineWidth', 2)
    plot(Data.head_free.x, Data.head_fixed.y, '-', 'MarkerSize', 10, 'LineWidth', 2)

end

%% Cluster data
function [y_med_sort] = cluster_data(I, K)
    [~, n_x] = size(I);
    y_med = nan(n_x,K);
    rng(1)
    for x = 1:n_x - 54
        C = I(:,x);
        y = find(C);
        if ~isempty(y)
            if K > 1
                g = kmeans(y, K); % sort into K groups
                y_group = cell(1,K);

                for k = 1:K
                   y_group{k} = y(g == k);
                   y_med(x,k) = median(y_group{k});
                end
            else % only one cluster
                y_med(x) = median(y);
            end
        end
    end
    y_med_sort = sort(y_med, 2); % [body , head , gaze]
end


%% Button function
function buttonPushed(src, event)
    disp('Extracting rectangle cooridnates')
    rect = findobj(gca, 'Type', 'images.roi.rectangle');
    n_rect = length(rect);
    
    % Get processed image reigons
    global H
    H = [];
    for n = 1:n_rect
        H.(rect(n).Label) = rect(n).Position;
    end
    
    close
end
