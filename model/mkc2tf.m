function [G] = mkc2tf(MM,KK,CC)
% mkc2tf: converts mass, damping, & stiffness matricies to transfer funciton matrix
% only for system of order n = 2
%
syms s

n = 2;
M = sym('m', n);
C = sym('c', n);
K = sym('k', n);
F = sym('f', [n 1]);
X = sym('x', [n 1]);
S0 = [s ; s].^0;
S1 = [s ; s].^1;
S2 = [s ; s].^2;

eom = collect(M*X.*S2 + C*X.*S1 + K*X.*S0 == F,X);

x2 = solve(eom(2), X(2));
EOM_1 = subs(eom(1), X(2), x2);
EOM_11 = subs(EOM_1, F(2), 0);
EOM_21 = subs(EOM_1, F(1), 0);
x11 = solve(EOM_11, X(1));
x21 = solve(EOM_21, X(1));

x1 = solve(eom(1), X(1));
EOM_2 = subs(eom(2), X(1), x1);
EOM_12 = subs(EOM_2, F(2), 0);
EOM_22 = subs(EOM_2, F(1), 0);
x12 = solve(EOM_12, X(2));
x22 = solve(EOM_22, X(2));

G(1,1) = collect(simplify(expand(x11/F(1))),s);
G(2,1) = collect(simplify(expand(x21/F(2))),s);
G(1,2) = collect(simplify(expand(x12/F(1))),s);
G(2,2) = collect(simplify(expand(x22/F(2))),s);

if nargin == 3
   G = subs(G, [M C K], [MM CC KK]);
end

end