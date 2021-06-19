function [C] = match_closedloop(p1, p2, h1, h2)
% match_closedloop: take plants and closed-loop models and designs
%                   controllers that shape plant into closed-loop model.
%
%  INPUTS:
%       p1  : plant #1
%       p2  : plant #2
%       h1  : plant #1
%       h2  : plant #2
%
%  OUTPUTS:
%       N   : numerator coeffcients
%

eq = controller_from_cl();

end