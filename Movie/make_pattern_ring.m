function [hpatch] = make_pattern_ring(pattern,pos,center,radius,thickness,color)
%% make_pattern_ring: draws the virtual world of the flight arena in the current axes
%   INPUT:
%       pattern   	: pattern structure
%                       Pats    : 4D panel map
%                       x_num   : number of x_positions
%                       gs_val  : grey scale level (0-3)
%       pos         : pattern position (x,y)
%       center      : center of pattern
%       radius    	: radius of pattern
%       thickness 	: thickness of pattern
%       color       : color of pattern
%   OUTPUT:
%       -
%

% pattern     = pattern_data.pattern;
% pos         = [50,5];
% center      = [round(FLY.yP/2) , round(FLY.xP/2)];
% radius      = floor(max([FLY.yP FLY.xP])/1.7);
% thickness   = 15;

if nargin<6
    color = 'g'; % default color
end

rout        = radius + thickness; % outer radius
res         = 360/pattern.x_num; % pattern resolution
sA          = deg2rad(res); % angle pixels subtend
theta_map   = deg2rad(res*(1:pattern.x_num)); % all pixel locations [rad]
Pats        = pattern.Pats(1,:,pos(1),pos(2)); % pixel map

% Draw pattern
hpatch = gobjects(length(theta_map),1);
hold on
for kk = 1:length(theta_map)
    xin = center(1) + radius*cos(theta_map(kk));
    xout = center(1) + rout*cos(theta_map(kk));
    xinN = center(1) + radius*cos(theta_map(kk) + sA);
    xoutN = center(1) + rout*cos(theta_map(kk) + sA);
    yin = center(2) + radius*sin(theta_map(kk));
    yout = center(2) + rout*sin(theta_map(kk));
    yinN = center(2) + radius*sin(theta_map(kk) + sA);
    youtN = center(2) + rout*sin(theta_map(kk) + sA);

	hpatch(kk) = patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin], color, ...
                'linestyle', 'none', 'FaceAlpha', Pats(kk)*(1/(2^(pattern.gs_val)-1)));
end

end
