function [] = run_ss_amp_norm(T, Fs, res, F, A, norm_vel, cent, root)
%% run_ss_amp_norm: makes normalized velocity ss functions
%
%   INPUT:
%       T           :   total time [s]
%       Fs          :   sampling frequency [Hz]
%       res         :   arena resolution
%       F           :   frequencies in signal
%       A           :   amplitude vector (+/-), can be salar for all same amplitude [°]
%       norm_vel   	:   peak velocity of each frequency component. "A" must be empty to set this. [°/s]
%       cent        :   center pixel on arena
%       showplot  	:   boolean (1= show time & freuency domain of signal)
%       root:       :   root directory to save position function file. Don't save if empty.
%
%
%   OUTPUT:
%       - 
%

clear ; close all ; clc

root        = 'C:\Users\BC\Box\Git\Fly-Multibody\Arena\functions';
F           = [];
T           = 20;
Fs          = 500;
A           = 3.75*[16 11 7 5 3 2 1];
norm_vel   	= 250;
res         = 3.75;
cent        = 45;
showplot    = true;

n_amp = length(A);
freq_all = nan(n_amp,1);
for a = 1:n_amp
    All = make_sos(T, Fs, res, F, A(a), norm_vel, cent, showplot, root);
    freq_all(a) = All.Freq;
end
disp(freq_all)

end