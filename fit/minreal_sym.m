function [sys] = minreal_sym(sys, var)
% minreal_sym:  normalized transfer function such that largest 
%            	denominator coefficient is equal to 1
%
%  INPUTS:
%       sys : symbolic transfer function. Must be linear.
%       var : symbolic variable. If not specified or empty, must only be one variable in sys.
%
%  OUTPUTS:
%       sys : minimum realization of system
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

% num = fliplr(coeffs(N, var));
den = fliplr(coeffs(D, var));

sys = (N/den(1)) / (D/den(1));

end