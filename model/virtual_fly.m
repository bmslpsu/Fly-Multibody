classdef virtual_fly
    %virtual_fly: construct a 2D animation of a fly
    %   Draws fly at given orientation
    %   Can be updated iteratively
    %
    
    properties (SetAccess = public, Hidden = false)
      	body            % body Ellipse object
        head            % head Ellipse object
        lwing           % left wing Arc object
     	rwing           % right wing Arc object
        centroid        %
    end
    
	properties (SetAccess = private, Hidden = false)
        wing_off        % wing offset ratio from centroid
        wing_off_L      % wing offset length from centroid
        head_centroid   % from semi-ellipse formula
    end
    
 	properties (SetAccess = public, Hidden = true)
    end

    methods
        function obj = virtual_fly(center, dim, head_color, body_color, alpha)
            %virtual_fly: Construct initial fly
            % 	center	: rotation point
            %  	dim  	: dimension structure
            %
            
            if nargin < 5
                alpha = 1;
                if nargin < 4
                    body_color = 'r';
                    if nargin < 3
                        head_color = 'b';
                        if nargin < 2
                            dim = 1;
                            if nargin < 1
                                center = [0 0];
                            end
                        end
                    end
                end
            end
            
            % Set the sizes
            if isstruct(dim)
                body_a = dim.a1;
                body_b = dim.b1;
                head_a = dim.a2;
                head_b = dim.b2;
                wing_L = 0.75*body_a;
                %rot_offset = -(body_a/2) - dim.d;
                rot_offset = dim.rot_off;
            else
                body_a = dim;
                body_b = 0.4*body_a;
                head_a = 0.4*body_a;
                head_b = 0.8*head_a;
                wing_L = 0.75*body_a;
                rot_offset = (body_a/2) - (body_a/7);
            end

            % Construct the fly
            default_body_ang = 0;
            default_head_ang = default_body_ang + 0;
            default_wing_ang = default_body_ang + 60;
            obj.body = Ellipse(center, body_a, body_b, 'ellipse', default_body_ang, rot_offset, body_color, alpha);
            obj.head = Ellipse(obj.body.top, head_a, head_b, 'semi', default_head_ang, 0, head_color, alpha);
            obj.rwing = Arc(obj.body.right, wing_L, default_wing_ang, 120, false, [0.3 0.1 0.7], 0.5);
            obj.lwing = Arc(obj.body.right, wing_L, 90 + default_wing_ang, 120, true, [0.3 0.1 0.7], 0.5);
            
            obj.wing_off = 0.3;
            obj.wing_off_L = obj.wing_off*obj.body.a/2;
            obj.head_centroid = obj.head.centroid + ( 2*obj.head.a / (3*pi) ) * ...
                [sind(obj.head.theta.a) cosd(obj.head.theta.a)];
            
            if ~nargin
                %obj = move(obj);
                obj = draw(obj,0,0);
            end
        end
        
        function obj = draw(obj, body_ang, head_ang, lwing_ang, rwing_ang, center)
            %DRAW Draw the shape
            %   Draw body & head
            %
            
            if nargin < 6
                center = obj.body.center;
                if nargin < 5
                    lwing_ang = obj.lwing.theta - obj.body.theta.a - 90;
                    if nargin < 4
                        rwing_ang = obj.rwing.theta - obj.body.theta.a;
                        if nargin < 3
                            head_ang = obj.head.theta.a - obj.body.theta.a;
                            if nargin < 2
                                body_ang = obj.body.theta.a;
                            end
                        end
                    end
                end
            end
            
            showpoints = false;
            
            % Update body
            obj.body.center = center;
            obj.body = update(obj.body, body_ang);
            
            % Update head
            obj.head.center = obj.body.top;
            obj.head = update(obj.head, body_ang + head_ang);
            
            % Update wings
            obj.rwing.center = obj.body.centroid + ...
                                    obj.wing_off_L*[sind(obj.body.theta.a) cosd(obj.body.theta.a)];
            obj.lwing.center = obj.rwing.center;
            
            % Draw fly
            obj.body  = draw(obj.body, body_ang, showpoints);
         	obj.rwing = draw(obj.rwing, rwing_ang - body_ang, showpoints);
         	obj.lwing = draw(obj.lwing, lwing_ang + body_ang, showpoints);
            obj.head  = draw(obj.head, body_ang + head_ang, showpoints);
            
            % Set colors
            set([obj.body.h.center], 'Marker', '.', 'MarkerEdgeColor', [0.6 0.6 0.6], ...
                'MarkerSize', 10*obj.body.a);
            uistack(obj.body.h.center, 'top')
            set([obj.body.h.centroid], 'Marker', '.', 'MarkerEdgeColor', [0.5 0 0], ...
                'MarkerSize', 25*obj.body.a);
            set([obj.body.h.long], 'Color', [0.7 0 0], 'LineWidth', 2)
            set([obj.head.h.long], 'Color', [0 0 0.7], 'LineWidth', 2)
            
            % Draw head centroid
            obj.head_centroid = obj.head.centroid + ( 2*obj.head.a / (3*pi) ) * ...
                [sind(obj.head.theta.a) cosd(obj.head.theta.a)];
            plot(obj.head_centroid(1), obj.head_centroid(2), '.', 'MarkerEdgeColor', [0 0 0.5], ...
                'MarkerSize', 35*obj.head.a);
        end
        
        function obj = move(obj,body_ang,head_ang,rwing,lwing,center)
            %MOVE Move the fly body & head
            %   Input the head and body angles
            %
            
            if nargin < 6
                Fs = 1000;
                t = (0:(1/Fs):10)';
                temp_body_ang = 1000*t + 50*sin(2*pi*15*t);
                %center = repmat(obj.body.center,length(t),1);
                center = 1*[(1:length(t))' , (1:length(t))'] ./ 2500;
                if nargin < 5
                    lwing = 60 + 15*sin(2*pi*30*t);
                   	if nargin < 4
                       	rwing = 60 + 15*cos(2*pi*30*t);
                        if nargin < 3
                            head_ang = 15*cos(2*pi*20*t);
                            if nargin < 2
                                body_ang = temp_body_ang;
                            end
                        end
                    end  
                end
            end
            
            assert(length(body_ang)==length(head_ang))
            N = length(body_ang);
            
            fig = figure;
            back_color = [204 255 255]./255;
            %set(fig,'Color',back_color)
            ax = subplot(1,1,1); hold on
            %set(ax, 'Color', back_color)
            axis image
            axis equal
            %axis off
            %move_range1 = obj.body.center + abs([obj.body.a , obj.body.a]);
            %move_range2 = obj.body.center - abs([obj.body.a , obj.body.a]);
            %axis([ sort([move_range1(1), move_range2(1)]) , sort([move_range1(2), move_range2(2)]) ])
            axis(max(center,[],'all')*[-1 1 -1 1])
            for n = 1:N
                obj = draw(obj, body_ang(n), head_ang(n), rwing(n), lwing(n), center(n,:));
                pause(1/Fs)
                cla
            end
            close
        end
        
    end
end

