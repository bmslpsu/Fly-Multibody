classdef system_stats
    % system_stats: calculates basic statistics for a matrix along a given
    % dimension ignoring nan's
    %   Calculates nanmean, nanmedian, nanstd
    
    properties (SetAccess = private, Hidden = false)
        mean
        median
        std
        circ_mean
        circ_median
        circ_std
    end
    
    methods
        function obj = system_stats(A,dim)
            % system_stats: Construct an instance of this class
            %   calculate mean, median, std
            
            if isnumeric(A)
                obj.mean        = nanmean(A,dim);
                obj.median      = nanmedian(A,dim);
                obj.std         = nanstd(A,[],dim);
                obj.circ_mean 	= circ_mean(A,[],dim);
                %obj.circ_median = circ_median(A,dim);
                obj.circ_std 	= circ_std(A,[],[],dim);
            else
                obj.mean        = A;
                obj.median      = A;
                obj.std         = A;
                obj.circ_mean 	= A;
                %obj.circ_median = A;
                obj.circ_std 	= A;
            end
        end
    end
end

