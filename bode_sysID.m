function [SYSTEM] = bode_sysID(time,ref,IOFv,debug,varargin)
%% bode_sysID: calculates weights of redundant actuation in frequency domain
%
%   INPUT:
%       time   	:   time vector
%       ref     :	reference input
%       IOFv   	:	input-output frequncies present in data (leave empty to automatically detect)
%
%   OUTPUT:
%       SYSTEM 	:   structure with system ID attributes
%
%   Usage: bode_sysID(time,ref,IOFv,input_1,input_2,...input_n)
%

nState  = length(varargin) + 1; % # of output states
nPoint  = length(time);         % # of data points in time domain
nFpoint = ceil(nPoint/2);       % # of data points in frequency domain
Fs      = 1/mean(diff(time));   % sampling rate [Hz]

State  = cellfun(@(x) x(:),varargin,'UniformOutput',false); % store inputs in cells as column vectors
State  = cat(2,State{:}); % store inputs in cells as column vectors
State(:,end+1)  = sum(State,2); % last state is the sum of all states 

% cmap        = jet(nOut);              % color map
fullColor   = [0.3 0.1 0.7];            % full state color
cmap        = [1 0 0;0 0 1;fullColor]; 	% color map
stimColor   = [0 1 0];                 	% stimulus color

% Convert reference to frequency domain
[refFv, refMag, refPhase, refFreq] = FFT(time,ref);
[IOFv, refIOMag, refIOPhase, IOidx] = getfreqpeaks(refFv, refMag, refPhase, IOFv, [], false);
refIOFreq = refFreq(IOidx);
nIO = length(IOFv);
nFreq = size(IOFv,1); % # of frequencies present
freqCmap = hsv(nFreq); % colormap across frequenncies

% Convert state outputs to frequency domain & calculate coherence
Mag     = nan(nFpoint,nState);
Phase   = nan(nFpoint,nState);
Freq 	= nan(nFpoint,nState);
Cohr   	= nan(nFpoint,nState);
IOMag  	= nan(nIO,nState);
IOPhase	= nan(nIO,nState);
IOFreq 	= nan(nIO,nState);
IOCohr 	= nan(nIO,nState);
for jj = 1:nState
    [Fv, Mag(:,jj), Phase(:,jj), Freq(:,jj)] = FFT(time,State(:,jj));
    [~, IOMag(:,jj),IOPhase(:,jj),IOidx] = getfreqpeaks(Fv, Mag(:,jj), Phase(:,jj), IOFv, [], false);
    IOFreq(:,jj) = Freq(IOidx,jj);
    [Cohr(:,jj),~] = mscohere(ref,State(:,jj),[],[],Fv,Fs);
    [~, IOCohr(:,jj),~,~] = getfreqpeaks(Fv, Cohr(:,jj), [], IOFv, [], false);
end

% Sort states by their overall means (to plot from largest to smallest)
IOmagMean = mean(IOMag,1);
[~,plotOrder] = sort(IOmagMean);
plotOrder = fliplr(plotOrder);

% Convert state outputs to frequency response functions
IOGain              = IOMag ./ refIOMag;
IOPhaseDiff         = IOPhase - refIOPhase;
IOFRF               = IOFreq ./ refIOFreq;
IOFRF_Gain          = abs(IOFRF);
IOFRF_PhaseDiff     = angle(IOFRF);

% IOPhaseDiff(IOPhaseDiff >  pi) = IOPhaseDiff(IOPhaseDiff >  pi) - 2*pi;
% IOPhaseDiff(IOPhaseDiff < -pi) = IOPhaseDiff(IOPhaseDiff < -pi) + 2*pi;
% 
% IOFRF_PhaseDiff(IOFRF_PhaseDiff >  pi) = IOFRF_PhaseDiff(IOFRF_PhaseDiff >  pi) - 2*pi;
% IOFRF_PhaseDiff(IOFRF_PhaseDiff < -4.8) = IOFRF_PhaseDiff(IOFRF_PhaseDiff < -4.8) + 2*pi;

% Calculate weights
IOWeight = IOGain(:,1:end-1) ./ IOGain(:,end);

SYSTEM.Mag              = Mag;
SYSTEM.Phase            = Phase;
SYSTEM.Cohr             = Cohr;
SYSTEM.IOMag            = IOMag;
SYSTEM.IOPhase      	= IOPhase;
SYSTEM.IOGain           = IOGain;
SYSTEM.IOPhaseDiff      = IOPhaseDiff;
SYSTEM.IOFRF            = IOFRF;
SYSTEM.IOFRF_Gain       = IOPhaseDiff;
SYSTEM.IOFRF_PhaseDiff  = IOPhaseDiff;
SYSTEM.IOCohr           = IOCohr;

if debug
    % Figure
    mrkSize = 10;
    fig(1) = figure; clf
    set(fig, 'Color', 'w','Units','inches','Position',[1 1 6 7])
    movegui(fig,'center')
    ax(1) = subplot(5,3,1:3); hold on
        href = plot(time, ref, 'Color', stimColor);
        hout = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(time, State(:,jj), 'Color', cmap(jj,:));
        end
        xlabel('Time (s)')
        ylabel('Magnitude (°)')
        leg = legend([href,hout],'Visual Motion','Body','Head','Gaze');
        leg.Box = 'off';
        leg.Orientation = 'horizontal';
        set(ax(1),'XLim',[0 time(end)])
        set([href hout],'LineWidth',1)

    ax(2) = subplot(5,3,4); hold on
        href = plot(refFv, refMag, 'Color', stimColor);
        hrefIO = plot(IOFv, refIOMag, '.-', 'Color', stimColor, 'MarkerSize', mrkSize);
        hout = gobjects(1,nState);
        houtIO = gobjects(1,nState);
        for jj = plotOrder
            hout(jj) = plot(Fv, Mag(:,jj), 'Color', cmap(jj,:));
            houtIO(jj) = plot(IOFv, IOMag(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylabel('Magnitude (°)')
        set([href hrefIO hout houtIO],'LineWidth',1)
        delete([hrefIO houtIO])

    ax(3) = subplot(5,3,7); hold on
        hrefIO = plot(IOFv, rad2deg(refIOPhase), '.-', 'Color', stimColor, 'MarkerSize', mrkSize);
        houtIO = gobjects(1,nState);
        for jj = 1:nState
            houtIO(jj) = plot(IOFv, rad2deg(IOPhase(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        set([hrefIO houtIO],'LineWidth',1)
        ylim(180*[-1 1])
        ylabel('Phase (°)')
        xlabel('Frequency (Hz)')

    ax(4) = subplot(5,3,5); hold on
        houtIO = gobjects(1,nState);
        for jj = 1:nState
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

    ax(5) = subplot(5,3,8); hold on
        houtIO = gobjects(1,nState);
        for jj = 1:nState
            houtIO(jj) = plot(IOFv, rad2deg(IOPhaseDiff(:,jj)), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        % ylim(180*[-1 1])
        ylim([-180 180])
        ylabel('Phase Difference(°)')
        xlabel('Frequency (Hz)')
        set([houtIO],'LineWidth',1)

    ax(6) = subplot(5,3,[6,9]); hold on
        houtIO = gobjects(1,nState-1);
        for jj = 1:nState-1
            houtIO(jj) = plot(IOFv, IOWeight(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylim([0 1.1])
        ylabel('Weight (°/°)')
        % xlabel('Frequency (Hz)')
        set([houtIO],'LineWidth',1)    

    ax(7) = subplot(5,3,[10,11,13,14]);
        houtIO = gobjects(1,nState);
        for jj = plotOrder
            houtIO(jj) = polarplot(IOFRF_PhaseDiff(:,jj), IOFRF_Gain(:,jj), '-', ...
                'Color', cmap(jj,:), 'MarkerSize', 1,'MarkerFaceColor','k','MarkerEdgeColor','k');
            hold on

            hfreq = gobjects(nFreq,1);
            for kk = 1:nFreq
                hfreq(kk) = polarplot(IOFRF_PhaseDiff(kk,jj), IOFRF_Gain(kk,jj), '.',...
                    'MarkerSize', 15,'MarkerFaceColor',freqCmap(kk,:),'MarkerEdgeColor',freqCmap(kk,:));
            end
        end

        ax(7) = gca;
    %     ax.RAxis.Label.String = 'Gain (°/°)';
    %     thetalabel('Phase Difference (°)')
        set(ax(7),'RLim',[0 1.1*maxGain])
        set([houtIO],'LineWidth',1)

        leg = legend(hfreq,string(IOFv));
        leg.Box = 'off';
        leg.Title.String = 'Frequency (Hz)';

    ax(8) = subplot(5,3,[12 15]); hold on
        hout = gobjects(1,nState);
        houtIO = gobjects(1,nState);
        for jj = 1:nState
            hout(jj) = plot(Fv, Cohr(:,jj), 'Color', cmap(jj,:));
            houtIO(jj) = plot(IOFv, IOCohr(:,jj), '.-', 'Color', cmap(jj,:), 'MarkerSize', mrkSize);
        end
        ylabel('Coherence')
        set([houtIO],'LineWidth',1)
        ylim([0 1.1])
        xlabel('Frequency (Hz)')
        % set(ax(2),'XLim',[0 time(end)])

    set(ax([2:6,8]), 'XLim', [0.2 1.1*max(IOFv)])
    set(ax([3,5]), 'YTick', -360:90:360)
    linkaxes(ax([2:6,8]),'x')
    
    set(ax,'LineWidth',1.5)
end

end