function [] = run_ss_vel_norm(T, Fs, res, F, A, norm_vel, cent, root)
%% run_ss_vel_norm: makes normalized velocity ss functions
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

root        = 'C:\Users\BC\Box\Git\Fly-Multibody\arena\functions';
T           = 10;
Fs          = 400.63;
res         = 3.75;
cent        = 45;
Phase       = -90;
showplot    = true;

A           = 3.75*[16 11 7 5 3 2 1];
norm_vel   	= 250;
F           = [];

% All = make_sos(T, Fs, res, F, A, norm_vel, cent, showplot, []);
% disp('Freq')
% disp(All.Freq)
% disp('Mean Vel')
% disp(mean(abs(All.dX)))
% disp('Median Vel')
% disp(median(abs(All.dX)))

n_amp = length(A); 
freq_all = nan(n_amp,1);
med_vel_all = nan(n_amp,1);
mean_vel_all = nan(n_amp,1);
for a = 1:n_amp
    All = make_sos(T, Fs, res, F, A(a), norm_vel, cent, Phase, showplot, []);
    freq_all(a) = All.Freq;
    med_vel_all(a) = median(abs(All.dX));
    mean_vel_all(a) = mean(abs(All.dX));
%     pause
%     close all
end
disp(freq_all)

end