function [] = batch_flipvid(root,varname,method,checkframe)
%% batch_flipvid: register each frame of a video to the same reference
%
%   INPUT:
%       root     	:   root directory
%       varname     :   video variable name
%       method      :   direction to flip
%                           'ud'        - flip vertically
%                           'lr'        - flip horizontally
%                           'rot180'	- flip vertically & horizontally
%       checkframe	:   visualize frame to determine if we flip the video (boolean)
%
%   OUTPUT:
%       -
%

if nargin < 4
    checkframe = true; % default
end

allmeth = ["ud","lr","rot90","rot180"];
methodI = find(strcmp(method, allmeth),1);
assert(~isempty(methodI), '"%s" is not a valid method', method)

% Select files
[files, path] = uigetfile('*.mat','Select file to switch', root, 'MultiSelect','on');
files = string(files);

% Flip videos
isWritable = true;
for kk = 1:length(files)
    disp(kk)
    flipframe = true;
    fname = fullfile(path,files(kk));
    matObj = matfile(fname, 'Writable', isWritable); % load .mat object
    if checkframe
        fig = figure ; clf ; cla
        vid = matObj.(varname);
        frame = vid(:,:,1,1);
        imagesc(frame) ; hold on; axis off
        title(files(kk),'Interpreter', 'none')
        uicontrol('Style','pushbutton','String','Flip','Callback',@flipcheck,'Position', [10 150 60 30])
        uicontrol('Style','pushbutton','String','Skip','Callback',@skipcheck,'Position', [10 190 60 30])
        waitfor(fig)
    end
    
    if flipframe
        disp('Flipping video...')
        matObj.flipped = method;
        switch methodI
            case 1  % flip video up/down
                matObj.(varname) = flipud(matObj.(varname));
          	case 2  % flip video left/right
                matObj.(varname) = fliplr(matObj.(varname));
          	case 3  % rotate video 90°
                matObj.(varname) = rot90(matObj.(varname), 1);
          	case 4  % rotate video 180°
                matObj.(varname) = rot90(matObj.(varname), 2);
            otherwise
                error('Not possible')
        end
    end
end
disp('Done')

% Callback functions for flip & skip buttons
%-------------------------------------------
function flipcheck(~,~)
    flipframe = true;
    close(fig)
end

function skipcheck(~,~)
    flipframe = false;
    close(fig)
end

end