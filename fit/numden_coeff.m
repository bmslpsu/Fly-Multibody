function [N,D] = numden_coeff(sys, var)
% numden_coeff: get the coefficicents of numerator and denominator of 
%               symbolic transfer function
%
%  INPUTS:
%       sys : symbolic transfer function. Must be linear.
%       var : symbolic variable. If not specified or empty, must only be one variable in sys.
%
%  OUTPUTS:
%       N   : numerator coeffcients
%    	D   : denominator coeffcients
%

if nargin < 2 % infer symbolic variable
    var = symvar(sys);
    var = unique(var);
    if isempty(var) % purely numeric transfer function, no computations needed
        N = double(sys);
        D = 1;
        return
    elseif length(var) > 1 % make sure there is only one variable
       error('symbolic TF may contain up to one symbolic variable')    
    end
else
    % use passed variable    
end

sys = collect(sys, var);
[N,D] = numden(sys);

N = fliplr(coeffs(N, var));
D = fliplr(coeffs(D, var));

end