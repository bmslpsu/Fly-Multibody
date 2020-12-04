function [mu, r_std, gain_mu, phase_mu, gain_std, phase_std] = complex_std(r, img, norm, rmv_out)
%% complex_std: standard deviation in the complex plane
%
%   INPUT:
%       r           :   real parts
%       img         :   imaginary parts
%       norm        :   normalize to mean or median
%       rmv_out     :   remove outliers before computing stats (boolean)
%
%   OUTPUT:
%       mu          : mean or median value
%      	r_std       : distance from mu std
%      	gain_std  	: gain std
%      	phase_std  	: gain std
%

if nargin < 4
    rmv_out = [];
end

r = r(:);
img = img(:);
if size(r)~=size(img)
    error('Error: X & X vectors must be the same size')
else
    r = r(:);
    img = img(:);
end

% Gain & phase
cmplx = r + 1i*img;
gain = abs(cmplx);
phase = angle(cmplx);

% Normalized to center
switch norm
    case 'mean'
        mu = [mean(r) mean(img)];
    case 'median'
        mu = [median(r) median(img)];
    otherwise
        error('"norm" must be "mean" or "median"')            
end
nr = r - mu(1);
nimg = img - mu(2);

% Find vectors from mean to points
R = sqrt(nr.^(2) + nimg.^(2));

if ~isempty(rmv_out)
    [R_rmv_out,I_out] = rmoutliers(R, 'grubbs');
    if length(R) ~= length(R_rmv_out)
        n_out = sum(I_out);
        warning('%i outliers detected', n_out)
    end
else
   R_rmv_out = R;
   I_out = false(length(R), 1);
end

% Mean/median gain & phase with outliers removed
gain_rmv_out = gain(~I_out);
phase_rmv_out = phase(~I_out);
switch norm
    case 'mean'
        mu = [mean(r(~I_out)) mean(img(~I_out))];
        gain_mu = mean(gain_rmv_out);
        phase_mu = circ_mean(phase_rmv_out);
    case 'median'
        mu = [median(r(~I_out)) median(img(~I_out))];
       	gain_mu = median(gain_rmv_out);
        phase_mu = circ_median(phase_rmv_out);    otherwise
        error('"norm" must be "mean" or "median"')            
end

% Calculate standard deviation of vectors from mean
r_std = std(R_rmv_out);

% Calculate standard deviation of gain % phase separately
gain_std = std(gain_rmv_out);
[~,phase_std] = circ_std(phase_rmv_out);

end