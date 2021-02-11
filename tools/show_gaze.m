function [h,top] = show_gaze(center,L,span,ang,C)
%% show_gaze: draw inverted triangle (gaze of fly)
%
%   INPUT:
%       center      : center point [x,y]
%       L           : length of gaze
%       span      	: field of view
%       ang       	: angle of ellispe in degrees
%       C       	: patch color
%
%   OUTPUT:
%       h           : handles for ellipse figure data
%       top         : top point of gaze
%       bot         : bottom point of gaze

% Example Input %
% center = [0,0];
% L = 10;
% span = 10;
% ang = 0;
% C = [0.2 0.1 0.5];

% Points
top     = center + L*[sind(ang),cosd(ang)];
left    = center + L*[sind(ang-span/2),cosd(ang-span/2)];
right 	= center + L*[sind(ang+span/2),cosd(ang+span/2)];

% All [x,y] points defining triangle
border  = [center;left;top;right];
xx      = border(:,1);
yy      = border(:,2);
vect    = [center;top];
xv      = vect(:,1);
yv      = vect(:,2);

% Plot
h = {};
hold on
h{end+1} = patch(xx,yy,C); % patch
alpha(h{end},0.2)
% h{end+1} = plot(xx,yy,'-*','Color',C);
% h{end+1} = plot(xv,yv,'-','Color',C,'LineWidth',1);

end