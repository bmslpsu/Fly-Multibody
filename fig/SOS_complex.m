function [] = SOS_complex()
%% SOS_magnitude:
root = 'E:\DATA\Magno_Data\Multibody';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'DATA','ALL','FUNC','GRAND','FLY','D','I','U','N')

%% FRF
clearvars -except DATA ALL GRAND FLY FUNC D I U N root
clc

% pI = [1 2 6 8 9 11 13];
% T = ["ref2body", "ref2head", "ref2wing", "head2body", "head2wing", "wing2body", "left2right"];

T = ["ref2body", "ref2head","ref2gaze"];
pI = [1 2 3];
yL = [1 0.6 1];

% T = ["ref2body_fixed"];
% pI = [1];
% yL = 1.2;

n_plot = length(pI);
n_freq = length(GRAND.all(1).IOFv(:,1,1));
fI = 2:n_freq;
n_freq_plot = length(fI);
cc = hsv(n_freq);

fig = figure (4) ; clf
% set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 9 3*n_plot])
mrksz = 6;

movegui(fig, 'center')
clear ax h hh
ax = gobjects(N{1,3},n_plot);
h = gobjects(N{1,3},n_plot, n_freq_plot);
leg = gobjects(N{1,3},n_plot);
IOFv = nan(n_freq_plot,N{1,3});


Gain = cell(n_plot,1);
Phase = cell(n_plot,1);
for v = 1:N{1,3}
    IOFv(:,v) = GRAND.fly_stats(v).mean.IOFv.mean(fI);
    all_complex = GRAND.all(v).IOFRF;
    G = abs(all_complex);
    P = angle(all_complex);
    R = real(all_complex);
    I = imag(all_complex);
    grand = GRAND.all_trial(v).IOFRF.median;
    %grand_std = GRAND.all_trial(v).IOFRF.std;
    for n = 1:n_plot
        subI = v + (n-1)*N{1,3};
        g = squeeze(G(:,pI(n),:));
        p = squeeze(P(:,pI(n),:));
        r = squeeze(R(:,pI(n),:));
        im = squeeze(I(:,pI(n),:));
        %ax(v,n) = subplot(N{1,3},n_plot,subI);
        ax(v,n) = subplot(n_plot,N{1,3},subI);
        [ax_cm,h_cm] = ComplexAxes(0:0.2:yL(n));
        set(ax_cm, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
        for f = 1:length(fI)
            plot(r(fI(f),:), im(fI(f),:), '.', 'MarkerSize', mrksz, ...
                    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 0.9*cc(f,:));
                        
            n_std = 2;
            [mu, r_std, gain_mu, phase_mu, gain_std, phase_std] = ...
                complex_std(r(fI(f),:), im(fI(f),:), 'median');
            
            h(v,n,f) = plot(mu(1), mu(2), '.', 'MarkerSize', 15, ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', cc(f,:));
            
            Gain{n}(f,v) = gain_std;
            Phase{n}(f,v) = phase_std;
            
            x = [0 mu(1)];
            y = [0 mu(2)];
            theta = rad2deg(atan2(x(2), y(2)));
            b = 2*pi*gain_mu*(rad2deg(phase_std) / 360);
            
            elps = Ellipse(mu, n_std*gain_std, n_std*b, ...
                'ellipse', theta, 0, cc(f,:), 0.4);
%             elps = Ellipse(mu, n_std*2*r_std, n_std*2*r_std, ...
%                 'ellipse', theta, 0, cc(f,:), 0.4);
            
            elps = draw(elps);
            set(elps.h.patch, 'EdgeColor', 'k')
            set(elps.h.centroid, 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
            hh.patch(v,n,f) = elps.h.patch;
            hh.med(v,n,f) = elps.h.centroid;
            
            %plot(x, y, '.k-', 'LineWidth', 1, 'MarkerSize', 10)
            %plot(mu(1), mu(2), '.', 'MarkerSize', 5, ...
                %'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
            
        end
        if n==1
            leg(v,n) = legend(squeeze(h(v,n,:)), string(IOFv(:,v)), 'Box', 'off', 'Location', 'east');
        end
        title(T(n))
        uistack(squeeze(hh.patch(v,n,:)), 'top')
        uistack(squeeze(h(v,n,:)), 'top')
    end
end

%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'Color', 'none', 'LineWidth', 1.5)


%% Save FRF data
% fname = 'Static_freq_mag';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end