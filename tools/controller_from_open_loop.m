function [C_new] = controller_from_open_loop(G, P_new)
% controller_from_open_loop: takes in the desired open-loop transfer
% function (G)and computes the required new controller (C_new) to shape a
% new plant (P_new) into G.
%
%  INPUTS:
%       G           : desired open-loop system
%       P_new   	: new plant
%
%  OUTPUTS:
%       sys         : structure containing model transfer functions
%      	data       	: bode plot and power data
%       h           : graphics handles
%

if isa(G,'sym')
    C_new = G / P_new;
else
    G_linear = tf([G.numerator{1}, G.denominator{1}]);
    C_new = G_linear / P_new;
end

end