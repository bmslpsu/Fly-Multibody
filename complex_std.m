function [mu, r_std, gain_mu, phase_mu, gain_std, phase_std] = complex_std(r, img, norm)
%% complex_std: standard deviation in the complex plane
%
%   INPUT:
%       r           :   real parts
%       img         :   imaginary parts
%       norm        :   normalize to mean or median
%
%   OUTPUT:
%       mu          : mean or median value
%      	r_std       : distance from mu std
%      	gain_std  	: gain std
%      	phase_std  	: gain std
%

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

% Normalized to mean
switch norm
    case 'mean'
        mu = [mean(r) mean(img)];
        gain_mu = circ_mean(gain);
        phase_mu = circ_mean(phase);
    case 'median'
        mu = [median(r) median(img)];
        gain_mu = circ_median(gain);
        phase_mu = circ_median(phase);
    otherwise
        error('"norm" must be "mean" or "median"')            
end
nr = r - mu(1);
nimg = img - mu(2);

% Find vectors from mean to points
R = sqrt(nr.^(2) + nimg.^(2));

% Calculate standard deviation of vectors from mean
r_std = std(R);

% Calculate standard deviation of gain % phase separately
gain_std = std(gain);
[~,phase_std] = circ_std(phase);

end