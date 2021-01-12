classdef singal_attributes
    % singal_attributes: calculates signal attributes
    %   
    
    properties (SetAccess = public, Hidden = false)
        Fs
        Ts
        time
        position
        velocity
        acceleration
        Fc
        position_lpf
        velocity_lpf
        acceleration_lpf
        
        detrend
        trend
        
        Fv
        mag
        phase
        freq
    end
    
    methods
        function obj = singal_attributes(position, time, Fc, n_detrend)
            % singal_attributes: Construct an instance of this class
            %   
            
            if nargin < 4
                n_detrend = [];
                if nargin < 3
                    Fc = [];
                end
            end
            
            obj.time = time;
            obj.Ts = mean(diff(obj.time));
            obj.Fs = 1 / obj.Ts;
            obj.position = position;
            obj.velocity = central_diff(obj.position, obj.Ts );
            obj.acceleration = central_diff(obj.position, obj.Ts );
            obj.Fc = Fc;
            
            if ~isempty(obj.Fc)
                [b,a] = butter(3, obj.Fc / (obj.Fs/2));
                obj.position_lpf = filtfilt(b, a, obj.position);
                obj.velocity_lpf = filtfilt(b, a, obj.velocity);
                obj.acceleration_lpf = filtfilt(b, a, obj.acceleration);
            end
            
            if ~isempty(n_detrend)
                obj.detrend = detrend(obj.position, n_detrend);
                obj.trend = obj.position - obj.detrend;
            end
            
            [obj.Fv, obj.mag.position, obj.phase.position, obj.freq.position] = FFT(obj.time, obj.position);
            [~, obj.mag.velocity, obj.phase.velocity, obj.freq.velocity] = FFT(obj.time, obj.velocity);
            [~, obj.mag.acceleration, obj.phase.acceleration, obj.freq.acceleration] = FFT(obj.time, obj.acceleration);
            
        end
    end
end

