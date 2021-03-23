classdef Arc
    %ARC construct a 2D ellipse shape using the patch function
    %   Draws shape at given orientation
    %   Can be updated iteratively
    %
    
    properties (SetAccess = public, Hidden = false)
      	center          % rotation point
        L               % dimension length 1 (long axis)
        theta           % orientation angle (°) CCW+ from right horizontal
      	offset          % rotation agle from orientation angle (°)
        lower           %
        flip            % flip the axis or not
        color           % patch face color
        facealpha       % pacth face opacity
     	h               % handles to graphics objects
    end 
    
	properties (SetAccess = private, Hidden = false)
        tip             % top coordinate
        axis            % long and short axes points
        border          % border points
       	span            % arc angle span
     	theta_global    % global coordinate system angle (°)
    end
    
 	properties (SetAccess = public, Hidden = false)
    end

    methods
        function obj = Arc(center,L,theta,offset,flip,color,facealpha)
            %ARC Construct an instance of this class
            %   Construct initial shape
            %
            
             % Defaults
            if nargin < 7
                facealpha = 0.2;;
                if nargin < 6
                    color = [0 0 0];
                    if nargin < 5
                        flip = false;
                        if nargin < 4
                            offset = 20;
                            if nargin < 3
                                theta = 0;
                                if nargin < 2
                                    L = 1;          
                                    if nargin < 1
                                        center = [0 0];
                                    end
                                end                                
                            end
                        end
                    end
                end
            end
            
            obj.center      = center;       % center point
            obj.L           = L;            % center point
            obj.theta       = theta;      	% orientation angle (°)
            obj.offset     	= offset;      	% rotation angle from orientation angle (°)
            obj.color       = color;       	% patch color
            obj.facealpha   = facealpha;  	% patch face opacity
            obj.flip        = flip;         % flip 0 axis by 180 & change to CW+
            
            % obj = draw(obj,theta,false);
        end
        
        function obj = update(obj,theta)
            %UPDATE Calculate geometric properties
            %   Grab current properties and use them to compute the geometric properties
            %
            
            obj.theta = theta;
            
            if ~obj.flip
                obj.theta_global = obj.theta;
                obj.lower = -obj.offset + obj.theta_global; % rotation angle from axis (°)
            elseif obj.flip
                obj.theta_global = 360-(obj.theta + 180);
                obj.lower = obj.offset + obj.theta_global; % rotation angle from axis (°)
            end
            
         	obj.tip = obj.center + obj.L*[cosd(obj.theta_global) , sind(obj.theta_global)]; % tip of arc
            obj.axis = [obj.center ; obj.tip]; % arc major axis
            obj.span = linspace(obj.theta_global,obj.lower, 100)'; % wrap around from orientation angle to constant lower angle
            obj.border = obj.center + obj.L*[cosd(obj.span) , sind(obj.span)]; % border in center reference frame
            obj.border = [obj.border ; obj.center ; obj.tip];
        end
        
        function obj = draw(obj,theta,showpoints)
            %DRAW Draw the shape
            %   Draw patch,centroid, center, border, surrounding points, axes
            %
            if nargin < 3
                showpoints = false;
                if nargin < 2
                    theta = obj.theta;
                end
            end
         	obj = update(obj,theta);
            
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

