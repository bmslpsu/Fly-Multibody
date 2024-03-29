function [All] = make_sos(T, Fs, res, F, A, norm_vel, cent, Phase, showplot, root)
%% make_sos: makes sum-of-sine position function
%
%   INPUT:
%     	T           :   total time [s]
%     	Fs          :   sampling frequency [Hz]
%     	res         :   arena resolution [�/pixel]
%     	F           :   vector containg frequencies in signal [Hz]
%      	A           :   amplitude vector (+/-), can be scalar for all frequencies with same amplitude [�]
%    	norm_vel   	:   peak velocity of each frequency component. "A" must be empty to set this. [�/s]
%    	cent        :   center pixel on arena
%      	Phase      	:   [�] phase vector, randomized if empty
%     	showplot  	:   boolean (1= show time & frequency domain of signal)
%   	root:       :   root directory to save position function file. Don't save if empty.
%
%   OUTPUT:
%       All         :   structure containing singal attributes
%

if nargin < 10
    root = [];
end

% Check inputs
assert( (length(T) == 1) && (T > 0), '"T" must be a positive scalar')
assert(all(F > 0), 'Frequencies must be positive')
assert( ((cent == 0) || (cent/round(cent) == 1)) && (cent >=0) && (cent <= 360/res), ...
    '"cent" must be a positive integer between 0 & 360/res')

freq_res = 1 / T; % frequency resolution
tt = (0:1/Fs:T)';  % time vector [s]
F = F(:);
A = A(:);
if isempty(A) % no amplitude set
    N = length(F); % # of frequency components
elseif isempty(F)
    N = length(A); % # of frequency components
else
    assert(length(A)==length(F), 'ampltidue and frequency vectors must be equal length if ''norm_vel'' is not specified')
    N = length(F);
end

if isempty(norm_vel)
    norm_vel = false;
else
    assert(isempty(A) || isempty(F), '"A" or "F" must be empty if "norm_vel" is specified')
end

if ~norm_vel
    if N == 1 % single-sine case
       assert(length(A) == length(N))
    else
        if length(A) == 1 % same amplitude for all frequencies
         	A = repmat(A,[N,1]);
        else % amplitude specified manually at each frequency
            assert(length(A) == N, ...
                'Specify "A" as a scalr or vector of length euqal to # of frequency components')
            A = A(:);
        end
    end
else
    if isempty(A)
        A = norm_vel ./ (2*pi*F); % normalize amplitude of each frequency component for given peak speed
    elseif isempty(F)
        F = norm_vel ./ (2*pi*A); % normalize frequency of each frequency component for given peak speed & ampltiude
    end
end
F = freq_res * round(F / freq_res); % round frequencies to nearest resolution point

% Check for prime multiples & shift up conflicting frequencies
prime_check = false(size(F));
while ~all(prime_check)
    Prime = nan(N);
    for n = 1:N
        for f = 1:N
            if (F(n) >= F(f))
                %P(n,f) = mod(F(n), F(f));
                Prime(n,f) = F(n) / F(f);
                check = rem(Prime(n,f),1);
                if (check < freq_res) && (F(n) ~= F(f))
                    newF = F(n) + freq_res;
                    disp(['Warning: ' num2str(F(n)) ...
                       ' is a prime multiple of ' num2str(F(f)), ...
                       '  /  change A = ' num2str(A(n)/3.75) ...
                       '  /  changing ' num2str(F(n)) ' to ' num2str(newF)])
                    F(n) = newF;
                    prime_check = false(size(F));
                else
                    prime_check(n) = true;
                end
            end
        end
    end
end

% Set phase
if isempty(Phase)
    disp('Randomizing phase')
    Phase = randi(359,N,1); % randomize initial phase of each frequency component [�]
elseif length(Phase) == 1 % same phase for all frequency components [�]
    Phase = repmat(Phase, [N,1]);
end
Phase = deg2rad(Phase); % [rad]

% Generate SOS signal
X = zeros(length(tt),1);
Y = zeros(length(tt),N);
for n = 1:N % each frequency component
    Y(:,n) = A(n)*cos(2*pi*F(n)*tt + Phase(n)); 
    X = X + Y(:,n);
end
dX = central_diff(X, 1/Fs);
dX_max = max(abs(dX));
dX_mean = mean(abs(dX));

X_step = res * round(X / res);
dX_step = central_diff(X_step, 1/Fs);
func  = round( (X_step / res) + cent); % function in panel position
% func  = (X / res) + cent; % function in panel position

% [Fv, mag.X , phase.X] = FFT(tt, X);
% [~, mag.dX , phase.dX] = FFT(tt, dX);
% [~, mag.X_step , phase.X_step] = FFT(tt, X_step);
% [~, mag.dX_step , phase.dX_step] = FFT(tt, dX_step);

% Fs_int = fix(Fs);
f1 = 0;
f2 = fix(Fs)/2;
[~,Fv,mag.X,phase.X] = chirpz(X, Fs, f1, f2);
[~,~,mag.dX,phase.dX] = chirpz(dX, Fs, f1, f2);
[~,~,mag.X_step,phase.X_step] = chirpz(X_step, Fs, f1, f2);
[~,~,mag.dX_step,phase.dX_step] = chirpz(dX_step, Fs, f1, f2);

% Assign outputs
All.func = func;
All.Fs = Fs;
All.time = tt;
All.res = res;
All.norm_vel = norm_vel;
All.Freq = F;
All.Amp = A;
All.Phase = Phase;
All.X = X;
All.dX = dX;
All.X_step = X_step;
All.dX_step = dX_step;
All.Fv = Fv;
All.mag = mag;
All.phase = phase;
All.freq_res = freq_res;

if showplot
    fig = figure (313);
    set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [6 2 10 8])
    clear ax h
    ax(1) = subplot(4,1,1) ; cla ; hold on; ylabel('Position (�)')
        h_smooth(1) = plot(tt, X, 'k');
        h_step(1) = plot(tt, X_step, 'r');
        ylim(1.1*max(abs(X))*[-1 1])

    ax(2) = subplot(4,1,2) ; cla ; hold on; ylabel('Velocity (�/s)') ; xlabel('Time (s)')
        h_smooth(2) = plot(tt, dX, 'k');
        %h_step(2) = plot(tt, dX_step, 'r');
        ylim(1.1*dX_max*[-1 1])

    ax(3) = subplot(4,1,3) ; cla ; hold on; ylabel('Position (�)')
        h_smooth(3) = plot(Fv, mag.X, 'k');
        h_step(2) = plot(Fv, mag.X_step, 'r');
        h_true(1) = plot(F, A, 'g.', 'MarkerSize', 12);
        legend([h_smooth(1) h_step(1) h_true(1)], 'continuous', 'discrete (arena)', 'true')

    ax(4) = subplot(4,1,4) ; cla ; hold on; ylabel('Velocity (�/s)') ; xlabel('Frequency (Hz)')
        h_smooth(4) = plot(Fv, mag.dX, 'k');  
        h_step(3) = plot(Fv, mag.dX_step, 'r');
        h_true(2) = plot(F, 2*pi*A.*F, 'g.', 'MarkerSize', 12);

    set(h_smooth, 'LineWidth', 1)
    set(h_step, 'LineWidth', 0.75)
    set(ax, 'LineWidth', 1.5,  'FontSize', 10)
    set(ax(1:2), 'XLim', [-0.5 T])
    set(ax(3:4), 'XLim', [0.9*min(F) 1.1*max(F)])
    %set(ax(3:4), 'XLim', [0 1.1*max(F)])
    set(ax([1,3]), 'XTickLabels', [])
    linkaxes(ax(1:2), 'x')
    linkaxes(ax(3:4), 'x')
end

% Save fucntion
if ~isempty(root)
    % Name file
    str_freq = '';
    str_amp = '';
	A_round = 0.05*round(A/0.05);
    for n = 1:N
        str_amp = [str_amp  num2str(0.01*round(A(n)/0.01)) '_'];
        if n < N
            str_freq = [str_freq  num2str(F(n)) '_'];
        else
            str_freq = [str_freq  num2str(F(n))];
        end
    end
    str_freq = strtrim(str_freq);
    if N == 1
        func_type = "SS";
    else
        func_type = "SOS";
    end
    if norm_vel
        fname = sprintf(['position_function_%s_Fs_%1.2f_T_%1.0f_vel_' num2str(norm_vel) ...
            '_amp_' str_amp 'freq_' str_freq '.mat'], func_type, Fs, T);
    else
        fname = sprintf(['position_function_%s_Fs_%1.2f_T_%1.0f' ...
            '_amp_' str_amp 'freq_' str_freq '.mat'], func_type, Fs, T);
    end
    fname_all = ['ALL_' fname];
    all_dir = fullfile(root, 'All');
    mkdir(all_dir)
    save(fullfile(root, fname), 'func')
    save(fullfile(all_dir, fname_all), 'All')
end

end