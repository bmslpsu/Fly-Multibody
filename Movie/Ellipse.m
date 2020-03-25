classdef Ellipse
    %ELLIPSE construct a 2D ellipse shape using the patch function
    %   Draws shape at given orientation
    %   Can be updated iteratively
    %
    
    properties (Access = public, Hidden = false)
      	center          % rotation point
        a               % dimension length 1 (long axis)
        b               % dimension length 2
      	offset          % area centroid offset ratio from centroid joint [da,db] between 0-1
        theta           % orientation angle (°)
        color           % patch face color
        facealpha       % pacth face opacity
     	h               % handles to graphics objects
    end 
    
	properties (Access = private, Hidden = false)
       	centroid       	% centroid coordinate
        top             % top coordinate
        bottom          % bottom coordinate
    	left            % left cooridinate
        right           % right cooridinate
        eccentricity    % eccentricity (ovalness)
        axis            % long and short axes points
        border          % border points
        border_normal   % border for 0°
        border_shift
        R               % rotation matrix for given theta
        w               % angle in set reference frame (CW positive measured from top vertical)
    end
    
 	properties (Access = public, Hidden = true)
        span
    end

    methods
        function obj = Ellipse(center,a,b,offset,theta,color,facealpha)
            %ELLIPSE Construct an instance of this class
            %   Construct initial shape
            %
            
            % Input checking
            if nargin < 7
                facealpha = 0.5; % default
                if nargin < 6
                    color = 'k'; % default
                    if nargin < 5
                        theta = 0; % default
                        if nargin < 4
                            offset = [0.5 0.5]; % default
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
            
            if b > a % a must be the long axis
                error('"a" must be the long axis')
            end
            
           	% Assign input properties
            obj.center   	= center;    	% center rotation point
            obj.a           = a;         	% long axis length
            obj.b           = b;           	% short axis length
            obj.theta.a     = theta;      	% long  axis direction
            obj.offset      = offset;     	% rotation point offset from centroid (ratio)
            obj.color       = color;       	% patch color
            obj.facealpha   = facealpha;  	% patch face opacity      
            
            % Short axis
         	obj.theta.b  = obj.theta.a + 90; % short axis direction
            

            
            obj = draw(obj,true);
        end
        
        function obj = draw(obj,showpoints)
            %DRAW Draw the shape
            %   Grab current properties and use them to compute the graphic properties
            %   Then draw the shape
            %
            
            if ~showpoints
                pcolor = 'none';
            else
            	pcolor = 'k';
            end
            
        	% Eccentricity
            obj.eccentricity = sqrt((obj.a^(2) - obj.b^(2))) / obj.a;
          	obj.offset = obj.offset + 0.5;
            
            % Top, bottom, centroid, left, right points
            obj.top      = obj.center + obj.a*(1-obj.offset(1))*[sind(obj.theta.a) , cosd(obj.theta.a)];
            obj.bottom   = obj.center - obj.a*(obj.offset(1))*[sind(obj.theta.a) , cosd(obj.theta.a)];  
            obj.centroid = (obj.top + obj.bottom) / 2;
           	obj.left	 = obj.centroid - obj.b*(1-obj.offset(2))*[sind(obj.theta.b) , cosd(obj.theta.b)];
            obj.right    = obj.centroid + obj.b*(obj.offset(2))*[sind(obj.theta.b) , cosd(obj.theta.b)];
            
            % Long & short axes
            obj.axis.a = [obj.top  ; obj.centroid ; obj.bottom]; % long  direction axis (a)
            obj.axis.b = [obj.left ; obj.centroid ; obj.right]; % short direction axis (b)
            
            % Rotation properties
         	obj.w = 180 - obj.theta.a; % rotation angle in set reference frame
            obj.R = [cosd(obj.w) -sind(obj.w); sind(obj.w) cosd(obj.w)]; % rotation matrix
            
            % Border
            obj.span            = linspace(0,2*pi,100)';
            obj.border_normal   = [obj.b/2*cos(obj.span) , obj.a/2*sin(obj.span)];        
            obj.border_shift    = (obj.R*obj.border_normal')'; % rotated border
            obj.border          = obj.centroid + obj.border_shift; % border in center reference frame
            
            % Draw the shape
            hold on
            % axis image
            obj.h.patch     = patch(obj.border(:,1), obj.border(:,2), obj.color, .... % patch
                                'FaceAlpha', obj.facealpha);
            obj.h.border    = plot(obj.border(:,1), obj.border(:,2), '-k', 'LineWidth', 1); % border
           	obj.h.long      = plot(obj.axis.a(:,1), obj.axis.a(:,2), '-k', 'LineWidth', 1); % long axis
            obj.h.short    	= plot(obj.axis.b(:,1), obj.axis.b(:,2), '-k', 'LineWidth', 1); % long axis
            obj.h.centroid 	= plot(obj.centroid(1), obj.centroid(2), '.k', ... % centroid
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 20);
            obj.h.center 	= plot(obj.center(1), obj.center(2), '.k', ... % center
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'g', 'MarkerSize', 30);    
         	obj.h.top       = plot(obj.top(1), obj.top(2), '.r', ... % top
                             	'MarkerFaceColor', 'none', 'MarkerSize', 20);
          	obj.h.bottom	= plot(obj.bottom(1), obj.bottom(2), '.k', ... % bottom
                             	'MarkerFaceColor', 'none', 'MarkerSize', 20);
         	obj.h.left     	= plot(obj.left(1), obj.left(2), '.m', ... % left
                             	'MarkerFaceColor', 'none', 'MarkerSize', 20);
         	obj.h.right    	= plot(obj.right(1), obj.right(2), '.k', ... % right
                             	'MarkerFaceColor', 'none', 'MarkerSize', 20);
        end
    end
end

