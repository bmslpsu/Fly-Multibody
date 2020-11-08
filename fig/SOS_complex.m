function [] = SOS_complex()
%% SOS_magnitude:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','GRAND','FLY','D','I','U','N')

%% FRF
clearvars -except DATA ALL GRAND FLY D I U N root
clc

pI = [1 2];
T = ["ref2body", "ref2head", "ref2gaze"];
n_plot = length(pI);
n_freq = length(GRAND.all(1).IOFv(:,1,1));
cc = hsv(n_freq);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
movegui(fig, 'center')
clear ax h
ax = gobjects(N.vel,n_plot);
h = gobjects(N.vel,n_plot, n_freq);
IOFv = nan(n_freq,N.vel);
for v = 1:N.vel
    all_complex = GRAND.all(v).IOFRF;
    G = abs(all_complex);
    P = angle(all_complex);
    IOFv(:,v) = GRAND.fly_stats(v).mean.IOFv.mean;
    for n = 1:n_plot
        subI = n + (v-1)*n_plot;
        g = squeeze(G(:,pI(n),:));
        p = squeeze(P(:,pI(n),:));
        ax(v,n) = subplot(N.vel,n_plot,subI);
        for f = 1:n_freq
            polarplot(p(f,:), g(f,:), '.', 'MarkerSize', 8, ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 0.7*cc(f,:));
                
            h(v,n,f) = polarplot(circ_median(p(f,:),2), median(g(f,:)), '.', 'MarkerSize', 20, ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', cc(f,:));
            hold on
        end
        title(T(n))
    end
end
uistack(h(v,n,f), 'top')
%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
%set(ax, 'LineWidth', 1.5)
legend(squeeze(h(1,1,:)), string(IOFv(:,1)), 'Box', 'off')

%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end