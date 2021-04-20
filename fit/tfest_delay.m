function [best_sys, sys_list, fitpercent, delay_sort] = tfest_delay(data, np, nz, opt, delay)
% tfest_delay: fits tranfer functions with set of time delays and picks best fit.
% Follows same input form as "tfest" but with extra input for delays to test
%
%  INPUTS:
%          data : system ID data object (frd, etc.)
%            np : # of poles
%            nz : # of zeros
%           opt : tfest options object
%         delay : delay vector
%
%  OUTPUTS:
%      best_sys : estimated transfer functions with best fit percent
%      sys_list : estimated transfer function list from best to worst fit percent
%    fitpercent : fit percents of each transfer fuction
%

% Defaults
if nargin < 5
    delay = 0; % acts like tfest with 0 delay
    if nargin < 4
       opt = tfestOptions; % use default options 
    end
end

% Fit for each delay
delay = delay(:);
n_delay = length(delay);
sys = cell(n_delay,1);
fitpercent = nan(n_delay,1);
for n = 1:n_delay
    sys{n} = tfest(data, np, nz, delay(n), opt);
    fitpercent(n) = sys{n}.Report.Fit.FitPercent;
end

% Sort fits and pull out best fit
[fitpercent, fit_order] = sort(fitpercent, 'descend');
delay_sort = delay(fit_order);
sys_list = sys(fit_order);
best_sys = sys_list{1};

end