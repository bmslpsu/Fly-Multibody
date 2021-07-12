function [cntrl, plant, tau, C, K, cutoff_freq] = get_controller(sys, I)
%% get_controller: seperates the plant and controller of a given transfer function 
% model based on a given a1 coefficient (inertia). Assumes plant and plant*controller
% are both first order, although numerator can be any order.
%
%   INPUT:
%     	sys     	: 1st order transfer function
%       I           : "inertia" of plant, if nan then normalize plant to have gain of one
%
%   OUTPUT:
%       cntrl   	: controller transfer function
%       plant   	: plant transfer function
%       C           : "damping" of plant
%       tau      	: time constant of plant [s]
%       K           : numerator of plant
%       cutoff_freq : cuttoff frequency of plant [hz]
%

sys = tf(sys);
num = sys.Numerator{1};
den = sys.Denominator{1};

assert(length(den) == 2, 'only 1st order transfer function models can be specified')

den_norm = den / den(2);
num_norm = num / den(2);
% sys_norm = tf(num_norm, den_norm);

tau = den_norm(1);
cutoff_freq = (1 / tau) / (2*pi);

if ~isnan(I)
    C = I / tau;
else
    C = 1;
end
K = 1 / C;

plant = tf(K, [tau 1]);
cntrl = num_norm / K;
cntrl = tf(cntrl,1);

end