classdef Ellipse
    %ELLIPSE construct a 2D ellipse shape using the patch function
    %   Draws shape at given orientation
    %   Can be updated iteratively
    %
    
    properties (SetAccess = public, Hidden = false)
        shape           % ellipse or semi-ellipse
      	center          % rotation point
        a               % dimension length 1 (long axis)
        b               % dimension length 2
      	offset          % area centroid offset ratio from centroid joint [da,db] between 0-1
        theta           % orientation angle (°)
        color           % patch face color
        facealpha       % pacth face opacity
     	h               % handles to graphics objects
    end 
    
	properties (SetAccess = private, Hidden = false)
       	centroid       	% centroid coordinate
        top             % top coordinate
        bottom          % bottom coordinate
    	left            % left cooridinate
        right           % right cooridinate
        eccentricity    % eccentricity (ovalness)
        axis            % long and short axes points
        border          % border points
        border_normal   % border at 0°
        border_shift    % rotated border
        R               % rotation matrix for given theta
        w               % angle in set reference frame (+CW measured from top vertical)
    end
    
 	properties (SetAccess = public, Hidden = true)
        span
        semi
    end

    methods
        function obj = Ellipse(center,a,b,shape,theta,offset,color,facealpha)
            %ELLIPSE Construct an instance of this class
            %   Construct initial shape
            %
            
            % Input checking
            if nargin < 8
                facealpha = 0.5; % default
                if nargin < 7
                    color = 'k'; % default
                    if nargin < 6
                        offset = 0; % default
                        if nargin < 5
                            theta = 0; % default
                            if nargin < 4
                                shape = 'ellipse';
                                if nargin == 2
                                    b = a; % default
                                elseif nargin < 2
                                    a = 1; % default
                                    b = 1; % default
                                    if ~nargin
                                        center = [0 0]; % default
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % % Error checking
            %if b > a % a must be the long axis
            %    error('"a" must be the long axis')
            %end
            
            if ~any( strcmp(shape,["ellipse","semi"]) ) % must be ellipse or semi-ellipse
                error('Accepted shapes are "ellipse" and "semi"')
            elseif strcmp(shape,"ellipse")
                obj.semi = false;
            elseif strcmp(shape,"semi")
                obj.semi = true;
            end
            
           	% Assign input properties
            obj.shape       = shape;        % shape
            obj.center   	= center;    	% center rotation point
            obj.a           = a;         	% long axis length
            obj.b           = b;           	% short axis length
            obj.theta.a     = theta;      	% long  axis direction
            obj.offset      = offset;     	% rotation point offset from centroid (ratio)
            obj.color       = color;       	% patch color
            obj.facealpha   = facealpha;  	% patch face opacity
            
            obj = update(obj);
            % obj = draw(obj,obj.theta.a,true);
        end
        
        function obj = update(obj,theta)
            %UPDATE Calculate geometric properties
            %   Grab current properties and use them to compute the geometric properties
            %
            
            if nargin == 2
                % Set new rotation angle
                obj.theta.a = theta;
            end
            
         	% Short axis
         	obj.theta.b = obj.theta.a + 90; % short axis direction
            
        	% Eccentricity
            obj.eccentricity = sqrt((obj.a^(2) - obj.b^(2))) / obj.a;
            off = obj.offset + 0.5;
            
            % Top, bottom, centroid, left, right points
            obj.top = obj.center + obj.a*(1-off)*[sind(obj.theta.a) , cosd(obj.theta.a)];
            if ~obj.semi
                obj.bottom = obj.center - obj.a*(off)*[sind(obj.theta.a) , cosd(obj.theta.a)];
                obj.centroid = (obj.top + obj.bottom) / 2;
            else
                obj.bottom = obj.center - obj.a*(off-0.5)*[sind(obj.theta.a) , cosd(obj.theta.a)];
                obj.centroid = obj.bottom;
            end
           	obj.left	= obj.centroid - obj.b*(1-0.5)*[sind(obj.theta.b) , cosd(obj.theta.b)];
            obj.right  	= obj.centroid + obj.b*(0.5)*[sind(obj.theta.b) , cosd(obj.theta.b)];
            
            % Long & short axes
            obj.axis.a = [obj.top  ; obj.centroid ; obj.bottom]; % long direction axis (a)
            obj.axis.b = [obj.left ; obj.centroid ; obj.right]; % short direction axis (b)
            
            % Rotation properties
            
         	obj.w = 180 - obj.theta.a; % rotation angle in set reference frame
            obj.R = [cosd(obj.w) -sind(obj.w); sind(obj.w) cosd(obj.w)]; % rotation matrix
            
            % Border
            if ~obj.semi
                obj.span = linspace(2*pi,0,100)'; % wrap around 0-2pi
            else
                obj.span = linspace(2*pi,pi,100)'; % wrap around 0-pi
            end
            obj.border_normal   = [obj.b/2*cos(obj.span) , obj.a/2*sin(obj.span)]; % border at 0°
            obj.border_shift    = (obj.R*obj.border_normal')'; % rotated border
            obj.border          = obj.centroid + obj.border_shift; % border in center reference frame
        end
        
        function obj = draw(obj,theta,showpoints)
            %DRAW Draw the shape
            %   Draw patch,centroid, center, border, surrounding points, axes
            %
            if nargin < 3
                showpoints = false;
            end
            
         	if nargin >= 2
             	obj = update(obj,theta);
          	end
            
            if ~showpoints
                pcolor = 'none';
            else
            	pcolor = 'k';
            end
            
            % Draw the shape
            hold on
            %axis image
            obj.h.patch     = patch(obj.border(:,1), obj.border(:,2), obj.color, .... % patch
                                'FaceAlpha', obj.facealpha, 'EdgeColor', pcolor);
                            
            obj.h.border    = plot(obj.border(:,1), obj.border(:,2), '-', 'Color', pcolor, ... % border
                                                'LineWidth', 1);
            
           	obj.h.long      = plot(obj.axis.a(:,1), obj.axis.a(:,2), '-', 'Color', pcolor, ... % long axis
                                                'LineWidth', 1);
                                            
            obj.h.short    	= plot(obj.axis.b(:,1), obj.axis.b(:,2), '-', 'Color', pcolor, ... % short axis
                                                'LineWidth', 1);
                                            
            obj.h.centroid 	= plot(obj.centroid(1), obj.centroid(2), '.', ... % centroid
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
                            
            obj.h.center 	= plot(obj.center(1), obj.center(2), '.', ... % center
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 30); 
                            
         	obj.h.top       = plot(obj.top(1), obj.top(2), '.', ... % top
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
                            
          	obj.h.bottom	= plot(obj.bottom(1), obj.bottom(2), '.', ... % bottom
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
                            
         	obj.h.left     	= plot(obj.left(1), obj.left(2), '.m', ... % left
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
                            
         	obj.h.right    	= plot(obj.right(1), obj.right(2), '.k', ... % right
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
        end
    end
end

