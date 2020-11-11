function [] = run_sos_custom(T, Fs, res, F, A, norm_vel, cent, root)
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
%   OUTPUT:
%       - 
%

clear ; close all ; clc

root        = 'C:\Users\BC\Box\Git\Fly-Multibody\Arena\functions';
F           = [];
T           = 20;
Fs          = 50;
res         = 3.75;
cent        = 45;
showplot    = true;

F           = [1 9.6];
N           = 5;
F           = linspace(F(1),F(2),N)'; % linearly spaced frequency vector [Hz]
F           = 0.1*round(F/0.1); % round frequencies to prime harmonics [Hz]
A           = [9;5;3;2;1];

All = make_sos(T, Fs, res, F, A, [], cent, showplot, root);

disp('Freq')
disp(All.Freq)
disp('Mean Vel')
disp(mean(abs(All.dX)))
disp('Median Vel')
disp(median(abs(All.dX)))
% histogram(All.dX)

end