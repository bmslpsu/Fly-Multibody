classdef FlyModel
    %FLYMODEL construct a 2D animation of a fly
    %   Draws fly at given orientation
    %   Can be updated iteratively
    %
    
    properties (SetAccess = public, Hidden = false)
      	body            % body Ellipse object
        head            % head Ellipse object
    end 
    
	properties (SetAccess = private, Hidden = false)
    end
    
 	properties (SetAccess = public, Hidden = true)
    end

    methods
        function obj = FlyModel(center,full_size,head_color,body_color,alpha)
            %FLYMODEL Construct an instance of this class
            %   Construct initial fly
            %
            
            % Set the sizes
            body_a = 0.75*full_size;
            body_b = 0.25*full_size;
            head_a = 0.25*full_size;
            head_b = 1.00*head_a;
            pin_offset = 0.5 - 1/7;
            
            % Construct the fly
            obj.body = Ellipse(center, body_a, body_b, 'ellipse', 0, pin_offset, body_color, alpha);
            obj.head = Ellipse(obj.body.top, head_a, head_b, 'semi', 0, 0, head_color, alpha);
            % obj = draw(obj,0,0);
        end
        
        function obj = draw(obj,body_ang,head_ang,center)
            %DRAW Draw the shape
            %   Draw body & head
            %
            
            if nargin < 4
                center = obj.body.center;
            end
            
            showpoints = false;
            
            obj.body.center = center;
            obj.body = draw(obj.body, body_ang, showpoints);
            
            obj.head.center = obj.body.top;
            obj.head = draw(obj.head, body_ang + head_ang, showpoints);
            
            set([obj.body.h.center,obj.head.h.bottom],'MarkerEdgeColor','k','MarkerSize',obj.head.a/2);
            % obj.head.h.bottom.MarkerEdgeColor = obj.head.h.patch.FaceColor;
        end
        
        function obj = move(obj,body_ang,head_ang,center)
            %MOVE Move the fly body & head
            %   Input the head and body angles
            %
            
            if nargin == 1
                t = 0:0.001:1;
                body_ang = 180*sin(2*pi*5*t);
                head_ang = 30*sin(2*pi*20*t);
            end
                    
            if nargin < 4
               center = repmat(obj.body.center,length(t),1);
               %center = [(1:length(t))' , (1:length(t))'];
            end
            
            assert(length(body_ang)==length(head_ang))
            N = length(body_ang);
            
            fig = figure;
            back_color = [204 255 255]./255;
            set(fig,'Color',back_color)
            ax = subplot(1,1,1); hold on
            set(ax,'Color',back_color)
            axis image
            axis equal
            axis off
            move_range1 = obj.body.center + abs([obj.body.a , obj.body.a]);
            move_range2 = obj.body.center - abs([obj.body.a , obj.body.a]);
            axis([ sort([move_range1(1), move_range2(1)]) , sort([move_range1(2), move_range2(2)]) ])
            % axis(max(center,[],'all')*[-1 1 -1 1])
            for n = 1:N
                obj = draw(obj,body_ang(n),head_ang(n),center(n,:));
                pause(0.001)
                cla
            end
            close
        end
        
    end
end

