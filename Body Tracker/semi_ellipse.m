function [h,top,bot] = semi_ellipse(center,L,ratio,ecc,ang,C)
%% draw_ellipse: draw an ellipse on figure window
%   INPUT:
%       center      : center point [x,y]
%       L           : length of ellipse
%       ratio       : portion of ellipse length below center
%       ecc       	: eccentricity
%       ang       	: angle of ellispe in degrees
%       C       	: patch color
%   OUTPUT:
%       h           : handles for ellipse figure data
%       top         : top point of ellipse
%       bottom    	: bottom point of ellipse

% Example Input %
% center = [0,0];
% L = 1;
% ratio = 0.5;
% ecc = 0.5;
% ang = 45;
% C = [1 0 0];

% Top point
top(1) = center(1) + 2*L*(1-ratio)*sind(ang);
top(2) = center(2) + 2*L*(1-ratio)*cosd(ang);

% Bottom point
bot(1) = center(1) - 2*L*ratio*sind(ang);
bot(2) = center(2) - 2*L*ratio*cosd(ang);

% All [x,y] points defining ellipse
xx = [center(1) top(1)];
yy = [center(2) top(2)];

% Calculate ellipse border points
a = 1/2*sqrt((top(1)-bot(1))^2+(top(2)-bot(2))^2);
b = a*sqrt(1-ecc^2);
t = linspace(pi/2,(3/2)*pi);
X = a*cos(t);
Y = b*sin(t);
w = atan2(bot(2)-top(2),bot(1)-top(1));
x = (top(1)+bot(1))/2 + X*cos(w) - Y*sin(w);
y = (top(2)+bot(2))/2 + X*sin(w) + Y*cos(w);
x = [x x(1)];
y = [y y(1)];

% Plot
h = {};
hold on
h{end+1} = patch(x,y,C); % patch
alpha(h{end},1)
h{end+1} = plot(x,y,'k-','LineWidth',2); % vector
h{end+1} = plot(xx,yy,'-k','LineWidth',2); % border
h{end+1} = plot(center(1),center(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % center point
% h{end+1} = plot(top(1),top(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % top point
% h{end+1} = plot(bot(1),bot(2),'-ko','MarkerSize',5,'MarkerFaceColor','k'); % bottom point
end