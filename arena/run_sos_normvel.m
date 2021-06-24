function [] = run_sos_normvel()
%% run_sos_amp_norm: makes normalized velocity ss functions
clear ; close all ; clc

root        = 'C:\Users\BC\Box\Git\Fly-Multibody\Arena\functions';
A           = [];
T           = 20;
Fs          = 400.63;
res         = 3.75;
cent        = 45;
phase       = [];
showplot    = true;

F = [0.4 1 2.4 3.5 5 8.3 11.1];
norm_vel = 50;
All = make_sos(T, Fs, res, F, A, norm_vel, cent, phase, showplot, []);

disp('Freq')
disp(All.Freq)
disp('Mean Vel')
disp(mean(abs(All.dX)))
disp('Median Vel')
disp(median(abs(All.dX)))
% histogram(All.dX)

end