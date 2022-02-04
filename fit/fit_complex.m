function [model, h] = fit_complex(cmplx_rsp, freq, num, den, tau, data_method, fit_method, showplot)
% fit_complex: fits tranfer function to emperical frequency response data
%
%  INPUTS:
%     cmplx_rsp : complex response, frequencies along rows & repetitions along columns
%          freq : frequency vector (hz)
%           num : polynomial form of numerator. ex: [1 1 1]
%           den : polynomial form of denominator. ex: [1 1 1]
%           tau : time constant. if empty, 0, or false, then model will not include a time constant
%               if true, then fit the time constant. if real, positive, scalar, fit using this constant
%               value
%   data_method : "RI" = fit real & imaginary parts, "GP" = fit gain & phase
%    fit_method : "lsqcurvefit" or "lsqnonlin"
%      showplot : show figure boolean
%
%  OUTPUTS:
%      model : structure containing model parameters and fit TF model
%          h : structure containing figure grpahics objects
%

% Defaults
if nargin < 8
    showplot = true;
    if nargin < 7
        fit_method = 'lsqcurvefit';
        if nargin < 6
            data_method = 'RI';
            if nargin < 5
                tau = [];
            end
        end
    end
end

% Check inputs
num = num(:)';
den = den(:)';
assert( all(any(num == [1 0]',1)) && all(any(den == [1 0]',1)), ...
    '"num" & "den" must be vectors of 1''s & 0''s''' )

% Fitting data 
if iscell(cmplx_rsp) % gain & phase in cells
    model.data.gain = cmplx_rsp{1};
    model.data.phase = cmplx_rsp{2};
    Cn = model.data.gain.* (cosd(model.data.phase) + 1i*sind(model.data.phase));
else % compelx response
    Cn = cmplx_rsp;
    model.data.gain = abs(Cn);
    model.data.phase = rad2deg(angle(Cn));
end
[n_freq, n_rep] = size(Cn);

% Convert frequency (hz) to radians
model.freq = freq;
model.rad = 2*pi*model.freq;

% Method parameters
model.data_method = data_method;
model.fit_method = fit_method;

% Model order
model.order.num = length(num) - 1;
model.order.den = length(den) - 1;
model.n_param_num = sum(num);
model.n_param_den = sum(den) - 1;

% Set required symbolic variables and make them real
syms s a0 b0 omega mm
assume([a0 b0 omega mm], 'real')

% Contruct plant numerator (a) and denominator (b) polynimals and make them real
model.sym.num = sym('a',[1,model.order.den+1]);
model.sym.num = fliplr([a0 model.sym.num(1:end-2) 1]);
model.sym.den = sym('b',[1,model.order.num+1]);
model.sym.den = fliplr([b0 model.sym.den(1:end-1)]);
assume([model.sym.num model.sym.den], 'real')
% a_coeff = model.sym.num(logical(den));
% b_coeff = model.sym.den(logical(num));

% Construct plant, with time delay if sepcified
model.sym.G = sum((num.*model.sym.den.*s.^(model.order.num:-1:0))) / sum((den.*model.sym.num.*s.^(model.order.den:-1:0)));

if isempty(tau) % no time constant, pass
    tc = 0;
elseif islogical(tau)
    if tau % fit the time constant
        model.sym.G = model.sym.G * exp(mm*s); % add time constant variable to plant
        tc = 1; % so we know to get the time constant from the fit
    else % no time constant, pass
        tc = 0;
    end
else % fit the model with specified time constant
    assert( (length(tau)==1) && isreal(tau) && (tau >= 0), 'tau must be a real, positive scalar')
    model.sym.G = model.sym.G * exp(tau*s);
    tc = 2;
end

% Substitute complex frequency into model (omega)
model.sym.Gjw = subs(model.sym.G, s, 1j*omega);

% Get real & imaginary parts and gain & phase of model
R = real(model.sym.Gjw);
I = imag(model.sym.Gjw);
Gain = sqrt(R^(2) + I^(2));
Phase = atan2(I,R);

% Convert real, imaginary, gain, & phase to anonymous functions
model.G = matlabFunction(model.sym.Gjw);
model.R = matlabFunction(R);
model.I = matlabFunction(I);
model.Gain = matlabFunction(Gain);
model.Phase = matlabFunction(Phase);
n_param = nargin(model.G) - 1; % # of parameters to fit

% Construct new anonymous functions for use with nonlinear curve-fitting minimization function
% Make string to evalute new anonymous functions
x_str = "";
for x = 1:n_param
    x_str = x_str + "x(" + string(x) + "), ";
end

% What data properties to fit
switch data_method
    case 'RI' % fit to real & imaginary parts
        X = [real(Cn(:)) , imag(Cn(:))]; % prepare data for fit
        fcn_str = ...
            sprintf('fnc = @(x,w) [model.R(%s w) , model.I(%s w)];', x_str, x_str);
        objfcn_str = ...
            sprintf('objfcn = @(x) [model.R(%s model.rad) - real(Cn), model.I(%s model.rad) - imag(Cn)];', x_str, x_str);
 	case 'GP' % fit to gain & phase
        X = [abs(Cn(:)) , angle(Cn(:))]; % prepare data for fit
        fcn_str = ...
            sprintf('fnc = @(x,w) [model.Gain(%s w) , model.Phase(%s w)];', x_str, x_str);
        objfcn_str = ...
            sprintf('objfcn = @(x) [model.Gain(%s model.rad) - abs(Cn), (model.Phase(%s model.rad) - angle(Cn))];', x_str, x_str);
    otherwise
       error('method must be "RI" or "GP"') 
end
eval(fcn_str) % construct anonymous function
eval(objfcn_str) % construct anonymous function

% Fit with positive intitial conditions but try negative if it doesn't work
% lb = -inf*ones(n_param,1)';
lb = 0.001*ones(n_param,1)';
ub = inf*ones(n_param,1)';
model.fitpercent.combined = 0;
x0_sign = [-1 1];
k = 1;
while (model.fitpercent.combined < 0.05) && (k <= 2)
    if k > 1
       warning('Trying new negative initial conditions')
    end
    
    % Set initial conditions and lower & upper bounds
    x0 = 0.05*x0_sign(k)*ones(n_param,1)';

    % What fitting algorithm to apply
    switch fit_method
        case 'lsqcurvefit'
            model.fcn = fnc;
            opts = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
            % x = lsqcurvefit(func.RI, x0, Fv, [real(Cn) , imag(Cn)], lb, ub, opts)
            freq_temp = repmat(model.rad, [size(Cn,2), 1]);
            x = lsqcurvefit(fnc, x0, freq_temp, X, lb, ub, opts);
        case 'lsqnonlin'
            model.fcn = objfcn;
            opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective', 'Display', 'off');
            %[x, resnorm, residuals, exitflag, output] = lsqnonlin(objfcn, x0, lb, ub, opts);
            x = lsqnonlin(objfcn, x0, lb, ub, opts);
        case 'fmincon'
    %         opts = optimset('TolFun', 10^-5, 'Display', 'none', 'MaxIter', 2000, 'MaxFunEvals', 4500);
    %         [x, V_eval, exitFlag(i), outputIteration] = ...
    %             fmincon(objfcn, x0, [], [], [], [], lb, ub,[], opts);
        otherwise
            error('method must be "lsqcurvefit" or "lsqnonlin"')
    end

    % Get fit parameters
    model.den = fliplr(x(1:model.n_param_den));
    model.num = fliplr(x(model.n_param_den+1:(model.n_param_den+1) + model.n_param_num-1));
    model.den = model.den .* den(2:end);
    model.num = model.num .* num;
    if tc == 1 % if time constant was fit
        model.tau = x(end);
    elseif tc == 0 % no time constant
        model.tau = 0;
    elseif tc == 2
        model.tau = tau;
    end

    % Construct the fit transfer function model
	model.G = tf(model.num, [1 model.den], 'IODelay', model.tau);

    % Check stability
 	stab = isstable(model.G);
    if ~stab
        warning('fit model is unstable')
    end

    % Get the real and imaginary trajectory and gain & phase
    [model.real, model.imag] = nyquist(model.G);
    model.real = squeeze(model.real);
    model.imag = squeeze(model.imag);

    [model.gain, model.phase, model.fvrad] = bode(model.G);
    model.gain = squeeze(model.gain);
    model.phase = squeeze(model.phase);
    model.fv = model.fvrad / (2*pi);

    % Get the real and imaginary trajectory  and gain & phase at input frequencies
    [model.real_freq, model.imag_freq] = nyquist(model.G, model.rad); 
    model.real_freq = squeeze(model.real_freq);
    model.imag_freq = squeeze(model.imag_freq);

    [model.gain_freq, model.phase_freq] = bode(model.G, model.rad);
    model.gain_freq = squeeze(model.gain_freq);
    model.phase_freq = squeeze(model.phase_freq);

    % Compute goodness of fit metrics
    gain_freq_all = repmat(model.gain_freq, [n_rep 1]);
    phase_freq_all = repmat(model.phase_freq, [n_rep 1]);
    gain_nmse = goodnessOfFit(model.data.gain(:), gain_freq_all, 'NMSE');
    phase_nmse = goodnessOfFit(model.data.phase(:), phase_freq_all, 'NMSE');
    model.fitpercent.gain = gain_nmse;
    model.fitpercent.phase = phase_nmse;
    model.fitpercent.combined = mean([model.fitpercent.gain model.fitpercent.phase]);

    k = k + 1;
end

if showplot
    cc = hsv(n_freq);
    freq_g = repmat(model.freq, [n_rep 1]);
    fig = figure; clf
    set(fig, 'Color', 'w', 'Units', 'inches', 'Visible', 'on')
    fig.Position(3:4) = [4.5 8.3];
    movegui(fig, 'center')
    ax(1) = subplot(3,1,1); cla ; hold on ; box on; set(ax(1), 'DataAspectRatio', [1 1 1])
        title('Nyquist')
        xline(0, '--k')
        yline(0, '--k')
        plot(1, 0, 'g.', 'MarkerSize', 25)
        h.model(1) = plot(model.real, model.imag, 'k');
        h.model_freq(1) = plot(model.real_freq, model.imag_freq, '.k', 'MarkerSize', 15);
        %plot(real(Cn), imag(Cn), '.r', 'MarkerSize', 10)
        h.data(:,1) = gscatter(real(Cn(:)), imag(Cn(:)), freq_g, cc, '.', 10, false);
        plot(model.real(1), model.imag(2), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 25)
        xlabel('Real')
        ylabel('Imaginary')
    ax(2) = subplot(3,1,2); cla ; hold on ; box on
        h.model(2) = plot(model.fv, model.gain, 'k');
        h.model_freq(2) = plot(model.freq, model.gain_freq, '.k', 'MarkerSize', 15);
        %plot(model.freq, abs(Cn), '.r', 'MarkerSize', 10);
        h.data(:,2) = gscatter(freq_g, model.data.gain(:), freq_g, cc, '.', 10, false);
        xlabel('Frequency (Hz)')
        ylabel('Gain')
    ax(3) = subplot(3,1,3); cla ; hold on ; box on
        h.model(3) = plot(model.fv, model.phase, 'k');
        h.model_freq(3) = plot(model.freq, model.phase_freq, '.k', 'MarkerSize', 15);
        %plot(model.freq, rad2deg(angle(Cn)), '.r', 'MarkerSize', 10)
        h.data(:,3) = gscatter(freq_g, model.data.phase(:), freq_g, cc, '.', 10, false);
        xlabel('Frequency (Hz)')
        ylabel('Phase')

    set(ax, 'Color', 'none', 'LineWidth', 1)
    set(ax, 'XGrid', 'on', 'YGrid', 'on')
    set(ax(2:3), 'XScale', 'log')
    linkaxes(ax(2:3), 'x')

    h.fig = fig;
    h.ax = ax;
else
    h = [];
end

end