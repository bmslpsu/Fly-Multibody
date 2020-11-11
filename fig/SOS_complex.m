function [] = SOS_complex()
%% SOS_magnitude:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','GRAND','FLY','FUNC','D','I','U','N')

%% FRF
clearvars -except DATA ALL GRAND FLY FUNC D I U N root
clc

pI = [2];
T = ["ref2body", "ref2head", "ref2gaze"];
n_plot = length(pI);
n_freq = length(GRAND.all(1).IOFv(:,1,1));
cc = hsv(n_freq);

fig = figure (2) ; clf
% set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 12 6])

movegui(fig, 'center')
clear ax h hh
ax = gobjects(N{1,3},n_plot);
h = gobjects(N{1,3},n_plot, n_freq);
leg = gobjects(N{1,3},n_plot);
IOFv = nan(n_freq,N{1,3});
for v = 1:N{1,3}
    IOFv(:,v) = GRAND.fly_stats(v).mean.IOFv.mean;
    all_complex = GRAND.all(v).IOFRF;
    G = abs(all_complex);
    P = angle(all_complex);
    R = real(all_complex);
    I = imag(all_complex);
    grand = GRAND.all_trial(v).IOFRF.median;
    %grand_std = GRAND.all_trial(v).IOFRF.std;
    for n = 1:n_plot
        subI = n + (v-1)*n_plot;
        g = squeeze(G(:,pI(n),:));
        p = squeeze(P(:,pI(n),:));
        r = squeeze(R(:,pI(n),:));
        im = squeeze(I(:,pI(n),:));
        %ax(v,n) = subplot(N{1,3},n_plot,subI);
        ax(v,n) = subplot(1,3,subI);
        [ax_cm,h_cm] = ComplexAxes(0:0.2:0.6);
        set(ax_cm, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
        for f = 1:n_freq        
            plot(r(f,:), im(f,:), '.', 'MarkerSize', 8, ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 0.65*cc(f,:));

            grand_r = real(grand(f,pI(n)));
            grand_im = imag(grand(f,pI(n)));
            %grand_theta = rad2deg(angle(grand(f,pI(n))));
                        
            n_std = 2;
            [mu, ~, gain_mu, ~, gain_std, phase_std] = ...
                complex_std(r(f,:), im(f,:), 'median');
            x = [0 mu(1)];
            y = [0 mu(2)];
            theta = rad2deg(atan2(x(2), y(2)));
            b = 2*pi*gain_mu*(rad2deg(phase_std) / 360);
            
            elps = Ellipse(mu, n_std*gain_std, n_std*b, ...
                'ellipse', theta, 0, cc(f,:), 0.4);
            elps = draw(elps);
            hh.patch(v,n,f) = elps.h.patch;
            hh.med(v,n,f) = elps.h.centroid;

            h(v,n,f) = plot(mu(1), mu(2), '.', 'MarkerSize', 15, ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', cc(f,:));
            
            %plot(x, y, '.k-', 'LineWidth', 1, 'MarkerSize', 10)
            %plot(mu(1), mu(2), '.', 'MarkerSize', 10 ...
                %'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
            
        end
        leg(v,n) = legend(squeeze(h(v,n,:)), string(IOFv(:,v)), 'Box', 'off');
        title(T(n))
    end
end
uistack(h(v,n,f), 'top')
uistack(hh.patch(v,n,f), 'top')
%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
%set(ax, 'LineWidth', 1.5)


%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end