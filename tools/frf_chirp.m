function [SYSTEM] = frf_chirp(time,ref,bin_sz,debug,varargin)
%% frf_chirp: calculates frequency response of chirp IO data
%
%   INPUT:
%       time        :   time vector or sampling frequency
%       ref         :	reference input
%       bin_sz    	:	frequency bin size (# of bins of size 1/T)
%       varargin    :   outputs
%
%   OUTPUT:
%       SYSTEM 	:   structure with system ID attributes
%           Names           - names of inputs states
%           Time            - time vector
%           State           - state matrix (one state per column)
%           Fv              - full frequency vector from chirp-z trasnform [Hz]
%           IOFv            - frequency vector at frequencies where FRF is defined [Hz]
%           refName         - reference (input) name
%           refState        - reference state
%           refFreq         - reference complex chirp-z trasnform output
%           refMag          - reference chirp-z trasnform magnitude
%           refPhase        - reference chirp-z transform phase [rad]
%           Freq            - state complex chirp-z trasnform output
%           Mag             - state chirp-z trasnform magnitude
%         	Phase           - state chirp-z transform phase [rad]
%         	FRF             - reference to state (input-output) complex response
%         	Gain            - state chirp-z transform phase at IO frequencies
%         	PhaseDiff       - state chirp-z transform phase at IO frequencies [rad]
%         	Cohr            - state coherence
%         	FRF_error       - state chirp-z error (from 1 + 0i)
%
%   Usage: SYSTEM = frf(time ,ref, bin_sz, input_1, input_2, ... , input_n)
%

% Store inputs in cells as column vectors
State = cellfun(@(x) x(:),varargin,'UniformOutput',false);
State = cat(2,State{:});
[nPoint,nState] = size(State); % # of data points in time domain & number of input states
nFpoint = ceil(nPoint/2); % # of data points in frequency domain

% Sampling rate [Hz]
if length(time)==1 % directly given
    Fs = time;   
else % from time vector
    Fs = 1 / mean(diff(time));
end

% Last state is the sum of all states if more than 1 state
if size(State,2) > 1
    State(:,end+1)  = sum(State,2);
    nState = nState + 1;
    
    % Get state names from input variables
    Names = string(nan(1,nState));
    for jj = 1:nState-1
        Names(jj) = inputname(4 + jj);
    end
    Names(nState) = "All";
else
    % Get state names from input variables
    Names = string(nan(1,nState));
    for jj = 1:nState
        Names(jj) = inputname(4 + jj);
    end
end

% Convert reference to frequency domain
refName = string(inputname(2));
[refFreq, refFv, refMag, refPhase] = chirpz(ref,Fs,0,Fs/2);

% Convert state outputs to frequency domain & calculate coherence
Mag     = nan(nFpoint,nState);
Phase   = nan(nFpoint,nState);
Freq 	= nan(nFpoint,nState);
Cohr   	= nan(nFpoint,nState);
for jj = 1:nState
    %[Fv, Mag(:,jj), Phase(:,jj), Freq(:,jj)] = FFT(time,State(:,jj));
    [Freq(:,jj), Fv, Mag(:,jj), Phase(:,jj)] = chirpz(State(:,jj),Fs,0,Fs/2);
    [Cohr(:,jj),~] = mscohere(ref,State(:,jj),[],[],Fv,Fs);
end

% Convert state outputs to frequency response functions
FRF = Freq ./ refFreq;

% Bin if specified
if bin_sz
    n_fv = length(FRF); % total frequency points
    FRF_binned = FRF; % new vector to store binned data
    for n = 1:nState
        for k = 1:n_fv
            win = (k-bin_sz):(k+bin_sz);
            win(win < 1) = [];
            win(win > n_fv) = [];
            FRF_binned(k,n) = mean(FRF_binned(win,n),1);
        end
    end
    FRF = FRF_binned;
end

% Extract gain & phase
% Gain        	= Mag ./ refMag;
% PhaseDiff    	= Phase - refPhase;
Gain          	= abs(FRF);
PhaseDiff   	= angle(FRF);
FRF_error    	= abs((1 + 0i) - FRF);

% lim1 = deg2rad(120);
% lim2 = deg2rad(300);
% 
% IOPhaseDiff(IOPhaseDiff >  lim1) = IOPhaseDiff(IOPhaseDiff >  lim1) - 2*pi;
% IOPhaseDiff(IOPhaseDiff < -lim2) = IOPhaseDiff(IOPhaseDiff < -lim2) + 2*pi;
% 
% IOFRF_PhaseDiff(IOFRF_PhaseDiff >  lim1) = IOFRF_PhaseDiff(IOFRF_PhaseDiff >  lim1) - 2*pi;
% IOFRF_PhaseDiff(IOFRF_PhaseDiff < -lim2) = IOFRF_PhaseDiff(IOFRF_PhaseDiff < -lim2) + 2*pi;
% 
% IOFRF_TimeDiff = (IOFRF_PhaseDiff ./ (2*pi)) .* repmat(1 ./ IOFv, [1 nState]);

% Calculate weights
% IOWeight = IOGain(:,1:end-1) ./ IOGain(:,end);
% IOWeight = IOGain(:,1:end-1) ./ sum(IOGain(:,1:end-1),2);
% if nState > 2
%     IORatio = IOWeight(:,2) ./ IOWeight(:,1);
% else
%     IORatio = ones(nIO,1);
% end

% Make output structure
SYSTEM.Names            = Names;
SYSTEM.Time             = time;
SYSTEM.State            = State;
SYSTEM.Fv               = refFv;
SYSTEM.refName          = refName;
SYSTEM.refState         = ref;
SYSTEM.refFreq          = refFreq;
SYSTEM.refMag        	= refMag;
SYSTEM.refPhase     	= refPhase;
SYSTEM.Freq             = Freq;
SYSTEM.Mag              = Mag;
SYSTEM.Phase            = Phase;
SYSTEM.FRF              = FRF;
SYSTEM.Gain             = Gain;
SYSTEM.PhaseDiff        = PhaseDiff;
SYSTEM.FRF_error        = FRF_error;
SYSTEM.Cohr             = Cohr;

% Figure
if debug
    % Sort states by their overall means (to plot from largest to smallest)
    magMean = mean(Mag,1);
    [~,plotOrder] = sort(magMean);
    plotOrder = fliplr(plotOrder);
    
    % Define colors
    fullColor   = [0.3 0.1 0.7];            % full state color
    cmap        = [1 0 0;0 0 1;fullColor]; 	% color map
    stimColor   = [0 1 0];                 	% stimulus color

    if nState ~=3
        cmap = [jet(nState-1) ; fullColor]; % color map
    end
    
    legLabel = [refName , Names];
    
    mrkSize = 10;
    fig(1) = figure; clf
    set(fig, 'Color', 'w','Units','inches','Position',[1 1 6 7])
    movegui(fig,'center')
    
    % Time Domain
    ax(1) = subplot(5,3,1:3); hold on
        href = plot(time, ref, 'Color', stimColor);
        hout = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(time, State(:,jj), 'Color', cmap(jj,:));
        end
        xlabel('Time (s)')
        ylabel('Magnitude (°)')
        leg = legend([href,hout],legLabel);
        leg.Box = 'off';
        leg.Orientation = 'horizontal';
        set(ax(1),'XLim',[0 time(end)])
        set([href hout],'LineWidth',1)
    
	% Frequency Magnitude
    ax(2) = subplot(5,3,4); hold on
        href = plot(refFv, refMag, 'Color', stimColor);
        hout = gobjects(1,nState);
        for jj = plotOrder
            hout(jj) = plot(Fv, Mag(:,jj), 'Color', cmap(jj,:));
        end
        ylabel('Magnitude (°)')
        set([href hout ],'LineWidth',1)
    
	% Frequency Phase
    ax(3) = subplot(5,3,7); hold on
%         hrefIO = plot(IOFv, rad2deg(refIOPhase), '.-', 'Color', stimColor, 'MarkerSize', mrkSize);
%         houtIO = gobjects(1,nState);
%         for jj = 1:nState
%             houtIO(jj) = plot(IOFv, rad2deg(IOPhase(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
%         end
%         set([hrefIO houtIO],'LineWidth',1)
        ylim(180*[-1 1])
        ylabel('Phase (°)')
        xlabel('Frequency (Hz)')

  	% IO Gain
    ax(4) = subplot(5,3,5); hold on
        hout = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(Fv, Gain(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        maxGain = max(max(Gain));
        if maxGain<1
            ylim([0 1.0])
        else
            ylim([0 1.1*maxGain])
        end
        ylabel('Gain (°/°)')
        set([hout],'LineWidth',1)
    
	% Phase Difference
    ax(5) = subplot(5,3,8); hold on
        hout = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(Fv, rad2deg(PhaseDiff(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        % ylim(180*[-1 1])
        ylim([-180 180])
        ylabel('Phase Difference(°)')
        xlabel('Frequency (Hz)')
        set([hout],'LineWidth',1)
    
    % Weights
    ax(6) = subplot(5,3,[6,9]); hold on
%         hout = gobjects(1,nState-1);
%         for jj = 1:nState-1
%             hout(jj) = plot(IOFv, IOWeight(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
%         end
%         ylim([0 1.1])
%         ylabel('Weight (°/°)')
%         % xlabel('Frequency (Hz)')
%         set([hout],'LineWidth',1)    

	% Complex Gain
    ax(7) = subplot(5,3,[10,11,13,14]);
        hout = gobjects(1,nState);
        for jj = plotOrder
            hout(jj) = polarplot(PhaseDiff(:,jj), Gain(:,jj), '.', ...
                'Color', cmap(jj,:), 'MarkerSize', 1, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
            hold on
        end

        ax(7) = gca;
        ax(7).RAxis.Label.String = 'Gain (°/°)';
        ax(7).ThetaAxis.Label.String = 'Phase Difference (°)';
        set(ax(7),'RLim',[0 1.1*maxGain])
        set([hout],'LineWidth',1)

	% Coherence & IO Coherence
    ax(8) = subplot(5,3,[12 15]); hold on
        hout = gobjects(1,nState);
        hout = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(Fv, Cohr(:,jj), 'Color', cmap(jj,:));
        end
        ylabel('Coherence')
        ylim([0 1.1])
        xlabel('Frequency (Hz)')

    set(ax([2:6,8]), 'XLim', [0.1 6])
    set(ax([3,5]), 'YTick', -360:90:360)
    linkaxes(ax([2:6,8]),'x')
    
    set(ax, 'LineWidth', 1.5)
end

end