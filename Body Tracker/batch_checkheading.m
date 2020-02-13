function [] = batch_checkheading(root)
%% batch_checkheading: loads in body tracked data & brings up initial figure frame to check for heading match
%
%   INPUT:
%       root   	:   root directory
%
%   OUTPUT:
%       -
%

% root = 'H:\EXPERIMENTS\MAGNO\Experiment_SOS\tracked_body';

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

fig = figure (1) ; clf
set(fig, 'Color', 'w')
hold on ; axis image
for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    
    load(fullfile(PATH,FILES(file)),'initframe')
    
    cla
    title(FILES(file), 'interpreter',' none') ; hold on
    imshow(initframe)
    
    pause
end
close all
disp('ALL DONE')
end