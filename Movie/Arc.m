classdef Arc
    %ARC construct a 2D ellipse shape using the patch function
    %   Draws shape at given orientation
    %   Can be updated iteratively
    %
    
    properties (SetAccess = public, Hidden = false)
      	center          % rotation point
        L               % dimension length 1 (long axis)
        theta           % orientation angle (°)
      	offset          % rotation agle from orientation angle (°)
        lower           %
        color           % patch face color
        facealpha       % pacth face opacity
     	h               % handles to graphics objects
    end 
    
	properties (SetAccess = private, Hidden = false)
        tip             % top coordinate
        axis            % long and short axes points
        border          % border points
       	span            % arc angle span
    end
    
 	properties (SetAccess = public, Hidden = false)
    end

    methods
        function obj = Arc(center,L,theta,offset,color,facealpha)
            %ARC Construct an instance of this class
            %   Construct initial shape
            %
            
            obj.center      = center;       % center point
            obj.L           = L;            % center point
            obj.theta       = theta;      	% orientation angle (°)
            obj.offset     	= offset;      	% rotation agle from orientation angle (°)
            obj.color       = color;       	% patch color
            obj.facealpha   = facealpha;  	% patch face opacity
            
            % obj = draw(obj,theta,false);
        end
        
        function obj = update(obj,theta)
            %UPDATE Calculate geometric properties
            %   Grab current properties and use them to compute the geometric properties
            %

            obj.theta = theta;
          	obj.lower = (-obj.offset + obj.theta); % rotation agle from orientation angle (°)
          	obj.tip = obj.center + obj.L*[cosd(obj.theta) , sind(obj.theta)]; % tip of arc
            obj.axis = [obj.center ; obj.tip]; % arc major axis
            obj.span = linspace(obj.theta, obj.lower, 100)'; % wrap around from orientation angle to constant lower angle
           	obj.border = obj.center + obj.L*[cosd(obj.span) , sind(obj.span)]; % border in center reference frame
            obj.border = [obj.border ; obj.center ; obj.tip];
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
            obj.h.patch     = patch(obj.border(:,1), obj.border(:,2), obj.color, .... % patch
                                'FaceAlpha', obj.facealpha, 'EdgeColor', pcolor);
                            
            obj.h.border    = plot(obj.border(:,1), obj.border(:,2), '-', 'Color', pcolor, ... % border
                                                'LineWidth', 1);
                                            
            obj.h.axis      = plot(obj.axis(:,1), obj.axis(:,2), '-', 'Color', pcolor, ... % border
                                                'LineWidth', 1);
                       
            obj.h.center 	= plot(obj.center(1), obj.center(2), '.', ... % center
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 30); 
                            
         	obj.h.tip       = plot(obj.tip(1), obj.tip(2), '.', ... % top
                             	'MarkerFaceColor', 'none', 'MarkerEdgeColor', pcolor, 'MarkerSize', 20);
        end
    end
end

