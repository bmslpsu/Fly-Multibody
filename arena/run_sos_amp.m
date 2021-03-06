function [] = run_sos_amp(T, Fs, res, F, A, norm_vel, cent, root)
%% run_sos_amp_norm: makes normalized velocity ss functions
%
%   INPUT:
%       T           :   total time [s]
%       Fs          :   sampling frequency [Hz]
%       res         :   arena resolution
%       F           :   frequencies in signal
%       A           :   amplitude vector (+/-), can be salar for all same amplitude [�]
%       norm_vel   	:   peak velocity of each frequency component. "A" must be empty to set this. [�/s]
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
Fs          = 400.63;
res         = 3.75;
cent        = 45;
phase       = 0;
showplot    = true;

% A = 1*[1 1 1 1 1 1 1 1];
% F = linspace(0.45,15, length(A))';

% A = 2*[1 1 1 1 1 1 1 1];
% F = linspace(0.35,13.05, length(A))';

% A = 3*[1 1 1 1 1 1 1 1];
% F = linspace(0.25, 11.85, length(A))';

% A = 1*[1 1 1 1 1 1 1 1];
% F = linspace(2.55, 19, length(A))';

% A = 2*[1 1 1 1 1 1 1 1];
% F = linspace(2.45, 18, length(A))';

% A = 1*ones(7,1);
% F = linspace(2.55, 17.55, length(A))';
% 
% A = 2*ones(7,1);
% F = linspace(2.4, 16.5, length(A))';
% 
% A = 3*ones(7,1);
% F = linspace(2.35, 15, length(A))';

A = 1*ones(7,1);
F = linspace(2.55, 15, length(A))';

A = 2*ones(7,1);
F = linspace(2.15, 13.05, length(A))';

A = 3*ones(7,1);
F = linspace(1.9, 11.85, length(A))';

All = make_sos(T, Fs, res, F, A, [], cent, phase, showplot, []);

disp('Freq')
disp(All.Freq)
disp('Mean Vel')
disp(mean(abs(All.dX)))
disp('Median Vel')
disp(median(abs(All.dX)))
% histogram(All.dX)

end