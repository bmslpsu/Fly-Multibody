function [controller] = get_controller(clsys, plant)
%% get_controller: solves for controller from closed-loop system transfer function and plant
%
%   INPUTS:
%       clsys 	:   closed-loop system transfer function
%       plant 	:   plant transfer function
%
%   OUTPUTS:
%       C       :   controller transfer function
%
syms s
clsys = tf(clsys);
G_N = clsys.Numerator{1};
G_D = clsys.Denominator{1};
G_nz = length(G_N)-1;
G_np = length(G_D)-1;

P_N = plant.Numerator{1};
P_D = plant.Denominator{1};
P_nz = length(P_N)-1;
P_np = length(P_D)-1;

[G,~,~,b,a] = symtf(G_nz, G_np, 'b', 'a');
[P,~,~,d,c] = symtf(P_nz, P_np, 'd', 'c');

C = G / ((1 - G)*P);
C = collect(simplify(expand(C)));

C = vpa(subs(C, [b a d c], [G_N G_D P_N P_D]),4);
[C_N,C_D] = numden(C);
scl = coeffs(C_D,s);
scl = scl(end);
C = (C_N/scl) / (C_D/scl);
C = sym2tf(C)

end