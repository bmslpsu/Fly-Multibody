function [] = run_sos_amp_norm(T, Fs, res, F, A, norm_vel, cent, root)
%% run_sos_amp_norm: makes normalized velocity ss functions
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
Fs          = 50;
%A        	= 3.75*[0.49 0.65 0.77 1.01 1.44 2 3 5 7 13];

% A        	= logspace((log(1.8)/log(10)),(log(45)/log(10)), 9)';
% norm_vel 	= 103;

A           = logspace((log(2.5)/log(10)),(log(61)/log(10)), 9)';
norm_vel   	= 148;

% A           = logspace((log(1.2)/log(10)),(log(61)/log(10)), 9)';
% norm_vel   	= 62;

res         = 3.75;
cent        = 45;
showplot    = true;

n_amp = length(A);
All = make_sos(T, Fs, res, F, A, norm_vel, cent, showplot, root);

disp('Freq')
disp(All.Freq)
disp('Mean Vel')
disp(mean(abs(All.dX)))
disp('Median Vel')
disp(median(abs(All.dX)))
% histogram(All.dX)

end