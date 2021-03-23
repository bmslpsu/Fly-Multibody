function [SYSTEM] = frf_v2(time, IOFv, input, outputs, min_coherence, showplot)
%% frf: calculates frequency response
%
%   INPUT:
%       time            : time vector or sampling frequency
%       IOFv            : input-output frequncies present in data (leave empty to automatically detect)
%       input           : input signal
%       outputs         : outputs [y1, y2, ..., yn]
%       min_coherence   : return nan for response under this coherence value
%       showplot        : show debug plot
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
%           refIOFreq       - reference complex chirp-z trasnform output at IO frequencies
%           refIOMag        - reference chirp-z trasnform magnitude at IO frequencies
%           refIOPhase      - reference chirp-z transform phase at IO frequencies [rad]
%           Freq            - state complex chirp-z trasnform output
%           Mag             - state chirp-z trasnform magnitude
%         	Phase           - state chirp-z transform phase [rad]
%           IOFreq          - state complex chirp-z trasnform output at IO frequencies
%           IOMag           - state chirp-z trasnform magnitude at IO frequencies
%         	IOPhase         - state chirp-z transform phase at IO frequencies [rad]
%         	FRF             - reference to state (input-output) complex response
%         	Gain            - state chirp-z transform phase at IO frequencies
%         	PhaseDiff       - state chirp-z transform phase at IO frequencies [rad]
%         	IOFRF           - reference to state (input-output) complex response at IO frequencies
%         	IOGain          - state chirp-z transform phase at IO frequencies
%         	IOPhaseDiff  	- state chirp-z transform phase at IO frequencies [rad]
%         	IOFRF_TimeDiff	- state chirp-z transform time difference (calculated from phase) at IO frequencies [rad]
%         	IOFRF_error     - state chirp-z error (from 1 + 0i) at IO frequencies
%         	Cohr            - state coherence
%         	IOCohr         	- state coherence at IO frequencies
%         	IOWeight      	- ratio of gains of input states over the total sum of the gains at IO frequencies
%         	IORatio        	- state coherence at IO frequencies
%

if nargin < 6
    showplot = false;
    if nargin < 5
        min_coherence = [];
    end
end

ftol = 0.01; % tolerance for finding peaks in frequency domain
IOFv = IOFv(:); % input-output frequncies
n_IO = length(IOFv); % # of IO frequencies

% Store inputs in cells as column vectors
State = outputs;
[n_point,n_state] = size(State); % # of data points in time domain & number of input states
n_fpoint = ceil(n_point/2); % # of data points in frequency domain

% Sampling rate [Hz]
if length(time)==1 % directly given
    Fs = time;   
else % from time vector
    Fs = 1 / mean(diff(time));
end

% Last state is the sum of all states if more than 1 state
if size(State,2) > 1
    State(:,end+1)  = sum(State,2);
    n_state = n_state + 1;
end

% Convert reference to frequency domain
[refFreq, refFv, refMag, refPhase] = chirpz(input,Fs,0,Fs/2);
[~, refIOMag, refIOPhase, IOidx] = getfreqpeaks(refFv, refMag, refPhase, IOFv, ftol, false);
refIOFreq = refFreq(IOidx);

% Convert state outputs to frequency domain & calculate coherence
Mag     = nan(n_fpoint,n_state);
Phase   = nan(n_fpoint,n_state);
Freq 	= nan(n_fpoint,n_state);
Cohr   	= nan(n_fpoint,n_state);
IOMag  	= nan(n_IO,n_state);
IOPhase	= nan(n_IO,n_state);
IOFreq 	= nan(n_IO,n_state);
IOCohr 	= nan(n_IO,n_state);
for jj = 1:n_state
    [Freq(:,jj), Fv, Mag(:,jj), Phase(:,jj)] = chirpz(State(:,jj),Fs,0,Fs/2);
    [~, IOMag(:,jj),IOPhase(:,jj),IOidx] = getfreqpeaks(Fv, Mag(:,jj), Phase(:,jj), IOFv, ftol, false);
    IOFreq(:,jj) = Freq(IOidx,jj);
    [Cohr(:,jj),~] = mscohere(input,State(:,jj),[],[],Fv,Fs);
    [~, IOCohr(:,jj),~,~] = getfreqpeaks(Fv, Cohr(:,jj), [], IOFv, ftol, false);
end

% Convert state outputs to frequency response functions
IOFRF        	= IOFreq ./ refIOFreq;
% IOGain      	= IOMag ./ refIOMag;
% IOPhaseDiff 	= IOPhase - refIOPhase;
IOGain          = abs(IOFRF);
IOPhaseDiff     = angle(IOFRF);
IOFRF_error  	= abs((1 + 0i) - IOFRF);

lim1 = deg2rad(120);
lim2 = deg2rad(300);
IOPhaseDiff(IOPhaseDiff >  lim1) = IOPhaseDiff(IOPhaseDiff >  lim1) - 2*pi;
IOPhaseDiff(IOPhaseDiff < -lim2) = IOPhaseDiff(IOPhaseDiff < -lim2) + 2*pi;
IOFRF_TimeDiff = (IOPhaseDiff ./ (2*pi)) .* repmat(1 ./ IOFv, [1 n_state]);

% Remove low coherence responses
if ~isempty(min_coherence)
    cohrI = IOCohr < min_coherence;
    IOFRF(cohrI) = nan;
    IOGain(cohrI) = nan;
    IOPhaseDiff(cohrI) = nan;
    IOFRF_error(cohrI) = nan;
    IOCohr(cohrI) = nan;
end

% % Calculate weights
% % IOWeight = IOGain(:,1:end-1) ./ IOGain(:,end);
% IOWeight = IOGain(:,1:end-1) ./ sum(IOGain(:,1:end-1),2);
% if n_state > 2
%     IORatio = IOWeight(:,2) ./ IOWeight(:,1);
% else
%     IORatio = ones(n_IO,1);
% end

% Make output structure
SYSTEM.Time             = time;
SYSTEM.State            = State;

SYSTEM.Fv               = refFv;
SYSTEM.IOFv             = IOFv;

SYSTEM.refState         = input;
SYSTEM.refFreq          = refFreq;
SYSTEM.refMag        	= refMag;
SYSTEM.refPhase     	= refPhase;
SYSTEM.refIOFreq     	= refIOFreq;
SYSTEM.refIOMag       	= refIOMag;
SYSTEM.refIOPhase    	= refIOPhase;

SYSTEM.Freq             = Freq;
SYSTEM.Mag              = Mag;
SYSTEM.Phase            = Phase;
SYSTEM.IOFreq         	= IOFreq;
SYSTEM.IOMag            = IOMag;
SYSTEM.IOPhase      	= IOPhase;

SYSTEM.IOGain           = IOGain;
SYSTEM.IOPhaseDiff      = IOPhaseDiff;
SYSTEM.IOFRF            = IOFRF;
SYSTEM.IOFRF_TimeDiff   = IOFRF_TimeDiff;
SYSTEM.IOFRF_error      = IOFRF_error;

SYSTEM.Cohr             = Cohr;
SYSTEM.IOCohr           = IOCohr;

% SYSTEM.IOWeight      	= IOWeight;
% SYSTEM.IORatio          = IORatio;

% Figure
if showplot
    % Sort states by their overall means (to plot from largest to smallest)
    IOmagMean = mean(IOMag,1);
    [~,plotOrder] = sort(IOmagMean);
    plotOrder = fliplr(plotOrder);
    
    % Define colors
    freqCmap    = hsv(n_IO); % colormap across frequenncies
    fullColor   = [0.3 0.1 0.7];            % full state color
    cmap        = [1 0 0;0 0 1;fullColor]; 	% color map
    stimColor   = [0 1 0];                 	% stimulus color

    if n_state ~=3
        cmap = [jet(n_state-1) ; fullColor]; % color map
    end
    
    %legLabel = [refName , Names];
    
    mrkSize = 10;
    fig(1) = figure; clf
    set(fig, 'Color', 'w','Units','inches','Position',[1 1 6 7])
    movegui(fig,'center')
    
    % Time Domain
    ax(1) = subplot(5,3,1:3); hold on
        href = plot(time, input, 'Color', stimColor);
        hout = gobjects(1,n_state);
        for jj = 1:n_state
            hout(jj) = plot(time, State(:,jj), 'Color', cmap(jj,:));
        end
        xlabel('Time (s)')
        ylabel('Magnitude (°)')
        %leg = legend([href,hout],legLabel);
        leg.Box = 'off';
        leg.Orientation = 'horizontal';
        set(ax(1),'XLim',[0 time(end)])
        set([href hout],'LineWidth',1)
    
	% Frequency Magnitude
    ax(2) = subplot(5,3,4); hold on
        href = plot(refFv, refMag, 'Color', stimColor);
        hrefIO = plot(IOFv, refIOMag, '.-', 'Color', stimColor, 'MarkerSize', mrkSize);
        hout = gobjects(1,n_state);
        houtIO = gobjects(1,n_state);
        for jj = plotOrder
            hout(jj) = plot(Fv, Mag(:,jj), 'Color', cmap(jj,:));
            houtIO(jj) = plot(IOFv, IOMag(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylabel('Magnitude (°)')
        set([href hrefIO hout houtIO],'LineWidth',1)
        delete([hrefIO houtIO])
    
	% Frequency Phase
    ax(3) = subplot(5,3,7); hold on
        hrefIO = plot(IOFv, rad2deg(refIOPhase), '.-', 'Color', stimColor, 'MarkerSize', mrkSize);
        houtIO = gobjects(1,n_state);
        for jj = 1:n_state
            houtIO(jj) = plot(IOFv, rad2deg(IOPhase(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        set([hrefIO houtIO],'LineWidth',1)
        ylim(180*[-1 1])
        ylabel('Phase (°)')
        xlabel('Frequency (Hz)')

  	% IO Gain
    ax(4) = subplot(5,3,5); hold on
        houtIO = gobjects(1,n_state);
        for jj = 1:n_state
            houtIO(jj) = plot(IOFv, IOGain(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        maxGain = max(max(IOGain));
        if maxGain<1
            ylim([0 1.0])
        else
            ylim([0 1.1*maxGain])
        end
        ylabel('Gain (°/°)')
        set([houtIO],'LineWidth',1)
    
	% Phase Difference
    ax(5) = subplot(5,3,8); hold on
        houtIO = gobjects(1,n_state);
        for jj = 1:n_state
            houtIO(jj) = plot(IOFv, rad2deg(IOPhaseDiff(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        % ylim(180*[-1 1])
        ylim([-180 180])
        ylabel('Phase Difference(°)')
        xlabel('Frequency (Hz)')
        set([houtIO],'LineWidth',1)
    
    % Weights
    ax(6) = subplot(5,3,[6,9]); hold on
        houtIO = gobjects(1,n_state-1);
        for jj = 1:n_state-1
            houtIO(jj) = plot(IOFv, IOWeight(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylim([0 1.1])
        ylabel('Weight (°/°)')
        % xlabel('Frequency (Hz)')
        set([houtIO],'LineWidth',1)    

	% Complex Gain
    ax(7) = subplot(5,3,[10,11,13,14]);
        houtIO = gobjects(1,n_state);
        for jj = plotOrder
            houtIO(jj) = polarplot(IOPhaseDiff(:,jj), IOGain(:,jj), '-', ...
                'Color', cmap(jj,:), 'MarkerSize', 1, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
            hold on

            hfreq = gobjects(n_IO,1);
            for kk = 1:n_IO
                hfreq(kk) = polarplot(IOPhaseDiff(kk,jj), IOGain(kk,jj), '.',...
                    'MarkerSize', 10, 'MarkerFaceColor', freqCmap(kk,:), 'MarkerEdgeColor', freqCmap(kk,:));
            end
        end

        ax(7) = gca;
        ax(7).RAxis.Label.String = 'Gain (°/°)';
        ax(7).ThetaAxis.Label.String = 'Phase Difference (°)';
        set(ax(7),'RLim',[0 1.1*maxGain])
        set([houtIO],'LineWidth',1)

        leg = legend(hfreq,string(IOFv));
        leg.Box = 'off';
        leg.Title.String = 'Frequency (Hz)';

	% Coherence & IO Coherence
    ax(8) = subplot(5,3,[12 15]); hold on
        hout = gobjects(1,n_state);
        houtIO = gobjects(1,n_state);
        for jj = 1:n_state
            hout(jj) = plot(Fv, Cohr(:,jj), 'Color', cmap(jj,:));
            %houtIO(jj) = plot(IOFv, IOCohr(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylabel('Coherence')
        %set([houtIO],'LineWidth',1)
        ylim([0 1.1])
        xlabel('Frequency (Hz)')
        % set(ax(2),'XLim',[0 time(end)])

    set(ax([2:6,8]), 'XLim', [0.2 5 + max(IOFv)])
    set(ax([3,5]), 'YTick', -360:90:360)
    linkaxes(ax([2:6,8]),'x')
    
    set(ax, 'LineWidth', 1.5)
end

end