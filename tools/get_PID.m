function [K] = get_PID(clsys, plant_order, dc_num, control_type)
%% get_PID: get PID gains based on transfer function of set order
%
%   INPUT:
%       sys             : transfer function object
%       plant_order     : order of plant numerator & denominator. if 2x1 then [den_order num_order]
%       dc_num          : dc magnitude or b0 of numerator
%       control_type   	: character string of any combination of 'P' , 'I' and 'D' to set type of controller
%
%   OUTPUT:
%       K           : PID gain structure
%

syms s a_1 a_0 K_p K_d K_i a0 b0

plant_order = plant_order(:);
assert((length(plant_order) <=2) && (length(plant_order) >=1), ...
    'plant order must a 1x1 or 1x2 vector')

assert((length(control_type) <=3) && ((length(control_type) >=1)), ...
    'control type must be ''PID'' or its compoenents')

if length(plant_order) == 1
    den_order = 2;
    num_order = 2;
elseif length(plant_order) == 2
    den_order = plant_order(1);
    num_order = plant_order(2);
end

% Contruct plant polynimals
a = sym('a',[1,den_order+1]);
a = fliplr([a0 a(1:end-2) 1]);

b = sym('b',[1,num_order+1]);
b = fliplr([b0 b(1:end-1)]);

if den_order==0
    a = a(2);
end

if dc_num && ~isempty(dc_num)
    b = subs(b, b0, dc_num);
end

% Plant
P = sum((b.*s.^(num_order:-1:0))) / sum((a.*s.^(den_order:-1:0)));

% PID controller
switch control_type
    case 'PID'
        % pass
    case 'PI'
        K_d = 0;
    case 'PD'
        K_i = 0;
    case 'P'
        K_d = 0;
        K_i = 0;
    case 'I'
        K_d = 0;
        K_p = 0;
    case 'D'
        K_p = 0;
        K_i = 0;
    otherwise
        error('control type must be ''PID'' or its compoenents')
end
C = K_d*s + K_p + K_i/s;

% Closed-loop system from plant controller derivation
G = P*C / (1 + P*C);
G = collect(simplify(expand(G)), s);
[N,D] = numden(G);

% Closed-loop system from identified model
% clsys_tf = tf(clsys);

N_coeff = fliplr(coeffs(N, s));
D_coeff = fliplr(coeffs(D, s));
cl_N_coeff = clsys.Numerator;
cl_D_coeff = clsys.Denominator;

% eq = [N_coeff == cl_N_coeff , D_coeff == cl_D_coeff];
Gains = solve(N_coeff == cl_N_coeff, [K_d K_p K_i]);

% G_cl = tf2sym(clsys_tf);
% G_cl = vpa(G_cl,3);
% % G_cl = vpa(collect(simplify(expand(G_cl))),3);
% [N_cl,D_cl] = numden(G_cl);

% n = length(cl_N_coeff);
% d = length(cl_D_coeff);
% eq = vpa([sum((cl_N_coeff.*s.^(n-1:-1:0))) == N,; ...
%           sum((cl_D_coeff.*s.^(d-1:-1:0))) == D],3);


end