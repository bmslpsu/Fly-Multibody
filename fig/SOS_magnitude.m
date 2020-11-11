function [] = SOS_magnitude()
%% SOS_magnitude:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','GRAND','FLY','D','I','U','N')

%% FFT
clearvars -except DATA ALL GRAND FLY D I U N root
clc

Fv = DATA.reference{1}.Fv;

fname = 'position';
% fname = 'velocity';

MAG.ref = cell(N.fly, N{1,3});
MAG.body = cell(N.fly, N{1,3});
MAG.head = cell(N.fly, N{1,3});
MAG.dwba = cell(N.fly, N{1,3});
for n = 1:N.file
    vel = DATA{:,3}(n);
    fly = DATA.fly(n);
    MAG.ref{fly,vel}(:,end+1) = DATA.reference{n}.mag.(fname);
    MAG.body{fly,vel}(:,end+1) = DATA.body{n}.mag.(fname);
    MAG.head{fly,vel}(:,end+1) = DATA.head{n}.mag.(fname);
    MAG.dwba{fly,vel}(:,end+1) = DATA.dwba{n}.mag.(fname);
end

% MAG = [];
% for v = 1:N{1,3}
%     for f = 1:N.fly
%        MAG.ref{f,v} = cat(2,ALL{f,v}.refMag);
%        MAG.body{f,v} = cat(3,ALL{f,v}.Mag);
%        MAG.head{f,v} = cat(2,ALL{f,v}.refMag);
%        MAG.dwba{f,v} = cat(2,ALL{f,v}.refMag);
%     end
% end

MAG.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    MAG, 'UniformOutput', false);

fnames = string(fieldnames(MAG.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:N{1,3}
        MAG.fly_mean.(fnames(f)){v} = cat(2, MAG.fly_stats.(fnames(f))(:,v).mean);
        MAG_all.(fnames(f)){v} = cat(2, MAG.(fnames(f)){:,v});
    end
end

MAG.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    MAG.fly_mean, 'UniformOutput', false);


%%
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 10 8])
ax = gobjects(n_field, N{1,3});
pp = 1;
cc = hsv(n_field);
for f = 1:n_field
    for v = 1:N{1,3}
        ax(f,v) = subplot(n_field, N{1,3}, pp);
            %xx = MAG.fly_mean.(fnames(f)){v};
            all = MAG_all.(fnames(f)){v};
            grand_mean = MAG.vel_stats.(fnames(f))(v).mean;
            grand_std = MAG.vel_stats.(fnames(f))(v).std;
            plot(Fv, all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            [~] = PlotPatch(grand_mean, grand_std, Fv, 0, 1, cc(v,:), 'k', 0.3, 1);
        
        pp = pp + 1;
    end
end

set(ax , 'LineWidth', 1.5, 'FontSize', 11, 'Box', 'off', 'XLim', [0.5 16])
% set(ax(1,:), 'YLim', [-1 160])
set(ax(2,:), 'YLim', [-0.1 2])
% set(ax(3,:), 'YLim', [-1 30])

linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')

set(ax(1:n_field-1,:), 'XTickLabel', [])
set(ax(1:n_field-1,:), 'XColor', 'none')
set(ax(:,2:N{1,3}), 'YTickLabel', [])
set(ax(:,2:N{1,3}), 'YColor', 'none')
set(ax(3,2:N{1,3}), 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Reference')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Body')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Head')
YLabelHC = get(ax(4,1), 'YLabel');
set([YLabelHC], 'String', '\DeltaWBA')
XLabelHC = get(ax(n_field,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')


%% Save magnitude data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end