function [model, h] = fitBode(cmplx_rsp, freq, num, den, tau, fit_method, showplot)
% fitBode: fits tranfer function to emperical frequency response data
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
if nargin < 7
    showplot = true;
    if nargin < 6
        fit_method = 'lsqcurvefit';
        if nargin < 5
            tau = [];
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
model.data.mean = structfun(@(x) mean(x,2), model.data, 'UniformOutput', false);
[n_freq, n_rep] = size(Cn);

% Convert frequency (hz) to radians
model.freq = freq;
model.rad = 2*pi*model.freq;

% Method parameters
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
    model.sym.G = model.sym.G * exp(-tau*s);
    tc = 2;
end

% Substitute complex frequency into model (omega)
model.sym.Gjw = subs(model.sym.G, s, 1j*omega);

% Get real & imaginary parts and gain & phase of model
model.sym.Real = real(model.sym.Gjw);
model.sym.Imag = imag(model.sym.Gjw);
model.sym.Gain = sqrt(model.sym.Real^(2) + model.sym.Imag^(2));
model.sym.Phase = atan2(model.sym.Imag,model.sym.Real);

% Convert real, imaginary, gain, & phase to anonymous functions
model.G = matlabFunction(model.sym.Gjw);
model.Real = matlabFunction(model.sym.Real);
model.Imag = matlabFunction(model.sym.Imag);
model.Gain = matlabFunction(model.sym.Gain);
model.Phase = matlabFunction(model.sym.Phase);
n_param = nargin(model.G) - 1; % # of parameters to fit

% Cost function
meanFlag = true;
switch meanFlag
    case true
        cost_all = sym(nan(n_freq, n_rep));
        for w = 1:n_freq
            for n = 1:n_rep
                cost_sym = log10(model.sym.Gain ./ model.data.gain(w,n)).^(2) + ...
                    (model.sym.Phase - deg2rad(model.data.phase(w,n))).^(2);
                cost_all(w,n) = subs(cost_sym, omega, model.rad(w));
            end
        end
    case false
        cost_all = sym(nan(n_freq,1));
        for w = 1:n_freq
            cost_sym = log10(model.sym.Gain ./ model.data.mean.gain(w)).^(2) + ...
                (model.sym.Phase - deg2rad(model.data.mean.phase(w))).^(2);
            cost_all(w) = subs(cost_sym, omega, model.rad(w));
        end 
end
cost_all = vpa(sum(cost_all,'all'),3);
cost_fnc = matlabFunction(cost_all);

% Construct new anonymous functions for use with nonlinear curve-fitting minimization function
% Make string to evalute new anonymous functions
x_str = "";
for x = 1:n_param
    x_str = x_str + "x(" + string(x) + "), ";
end
x_str = char(x_str);
x_str = x_str(1:end-2);

objfcn_str = ...
    sprintf('objfcn = @(x) cost_fnc(%s);', x_str);
eval(objfcn_str) % construct anonymous function

% Fit with positive intitial conditions but try negative if it doesn't work
% lb = -inf*ones(n_param,1)';
lb = 0.05*ones(n_param,1)';
ub = inf*ones(n_param,1)';
if tc == 1
   lb(end) = eps;
   ub(end) = 0.1;
end
model.fitpercent.combined = 0;
x0_sign = [1 -1];
k = 1;
while (model.fitpercent.combined < 0.05) && (k <= 2)
    if k > 1
       warning('Trying new negative initial conditions')
    end
    
    % Set initial conditions and lower & upper bounds
    x0 = 0.05*x0_sign(k)*ones(n_param,1)';

    % What fitting algorithm to apply
    switch fit_method
        case 'lsqnonlin'
            %model.fcn = objfcn;
            opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective', 'Display', 'off');
            %[x, resnorm, residuals, exitflag, output] = lsqnonlin(objfcn, x0, lb, ub, opts);
            x = lsqnonlin(objfcn, x0, lb, ub, opts);
        case 'fmincon'
            opts = optimset('TolFun', 10^-5, 'Display', 'none', 'MaxIter', 2000, 'MaxFunEvals', 4500);
            [x, V_eval, exitFlag, outputIteration] = ...
                fmincon(objfcn, x0, [], [], [], [], lb, ub, [], opts);
        otherwise
            error('method must be "lsqcurvefit" or "lsqnonlin"')
    end

    % Get fit parameters    
	model.num = fliplr(x(model.n_param_den+1:(model.n_param_den+1) + model.n_param_num-1));
    num_temp = num;
    num_temp(logical(num)) = model.num;
    model.num = num_temp;   
    
	model.den = fliplr(x(1:model.n_param_den));
    den_temp = den(2:end);
    den_temp(logical(den(2:end))) = model.den;
    model.den = den_temp;
    
    if tc == 1 % if time constant was fit
        model.tau = x(end);
    elseif tc == 0 % no time constant
        model.tau = 0;
    elseif tc == 2
        model.tau = tau;
    end

    % Construct the fit transfer function model
	model.tf.G = tf(model.num, [1 model.den], 'IODelay', model.tau);

    % Check stability
 	stab = isstable(model.tf.G);
    if ~stab
        warning('fit model is unstable')
    end

    % Get the real and imaginary trajectory and gain & phase
    [model.real, model.imag] = nyquist(model.tf.G);
    model.real = squeeze(model.real);
    model.imag = squeeze(model.imag);

    [model.gain, model.phase, model.fvrad] = bode(model.tf.G);
    model.gain = squeeze(model.gain);
    model.phase = squeeze(model.phase);
    model.fv = model.fvrad / (2*pi);

    % Get the real and imaginary trajectory  and gain & phase at input frequencies
    [model.real_freq, model.imag_freq] = nyquist(model.tf.G, model.rad); 
    model.real_freq = squeeze(model.real_freq);
    model.imag_freq = squeeze(model.imag_freq);

    [model.gain_freq, model.phase_freq] = bode(model.tf.G, model.rad);
    model.gain_freq = squeeze(model.gain_freq);
    model.phase_freq = squeeze(model.phase_freq);

    % Compute goodness of fit metrics
    gain_freq_all = repmat(model.gain_freq, [n_rep 1]);
    phase_freq_all = repmat(model.phase_freq, [n_rep 1]);
    gain_nmse = goodnessOfFit(model.data.gain(:), gain_freq_all, 'NMSE');
    phase_nmse = goodnessOfFit(model.data.phase(:), phase_freq_all, 'NMSE');
    model.fitpercent.gain = 1 - gain_nmse;
    model.fitpercent.phase = 1 - phase_nmse;
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