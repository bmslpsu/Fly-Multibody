function [A,B] = mkc2ss(MM,KK,CC)
%
n = 2;
M = sym('m', n);
C = sym('c', n);
K = sym('k', n);
F = sym('f', [n 1]);
X = sym('x', [3*n 1]);

eom = M*X(5:6) + C*X(3:4) + K*X(1:2) == F;
x6 = solve(eom(2), X(6));
x5 = X(5) == solve(eom(1), X(5));
x5 = subs(x5, X(6), x6);
x5 = solve(x5, X(5));
x5 = simplify(expand(x5));
EOM1 = collect(x5, [X ; F]);
[C1,T1] = coeffs(EOM1, [X ; F]);

x5 = solve(eom(1), X(5));
x6 = X(6) == solve(eom(2), X(6));
x6 = subs(x6, X(5), x5);
x6 = solve(x6, X(6));
x6 = simplify(expand(x6));
EOM2 = collect(x6, [X ; F]);
[C2,T2] = coeffs(EOM2, [X ; F]);
A = [0 1 0 0 ; C1([3 1 4 2]) ; 0 0 0 1 ; C2([3 1 4 2])];
B = [0 0; C1(5:6) ; 0 0 ; C2(5:6)];

if nargin == 3
   A = subs(A, [M C K], [MM CC KK]);
   B = subs(B, [M C K], [MM CC KK]); 
end

end