function [] = SOS_complex_free_fixed()
%% SOS_complex_free_fixed:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL');

%% Body-free head vs body-fixed head
clearvars -except ALL FILE PATH
clc

ALL.HeadFree.FRF_data.ref2head.fly(1).complex  

T = ["ref2body", "ref2head","ref2head"];
pI = [1 2 3];
yL = [1 1 1];

n_plot = length(pI);
n_cond = N{1,3};
% n_cond = 2;
% n_freq = length(GRAND.all(1).IOFv(:,1,1));
% fI = 1:n_freq;
% n_freq_plot = length(fI);
% cc = hsv(n_freq);

fig = figure (1) ; clf
% set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 3*n_cond])
mrksz = 6;

movegui(fig, 'center')
clear ax h hh
ax = gobjects(n_cond, n_plot);
h = gobjects(n_cond, n_plot, 9);
leg = gobjects(1, n_cond);
for n = 1:n_plot
    for v = 1:n_cond
        IOFv = FRF_data.IOFv{v};
        R = real(FRF_data.(T(n)).fly(v).complex(:));
        I = imag(FRF_data.(T(n)).fly(v).complex(:));
        [n_freq, n_fly] = size(FRF_data.(T(n)).fly(v).complex);
        G = repmat(IOFv, [n_fly, 1]);
        cc = hsv(n_freq);
        
        subI = v + (n-1)*n_cond;
        ax(v,n) = subplot(n_cond,n_plot,subI);
        %ax(v,n) = subplot(n_plot,n_cond,subI);
        [ax_cm,h_cm] = ComplexAxes(0:0.5:yL(n));
        set(ax_cm, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
        
        h(v,n,:) = gscatter(R, I, G, cc, '.', 7, false);
        mu = FRF_data.(T(n)).grand_mean(v).complex;
        plot(real(mu), imag(mu), '-k.', 'MarkerSize', 12, 'MarkerFaceColor', 'none');
            
        for f = 1:n_freq          
            x = [0 real(mu(f))];
            y = [0 real(mu(f))];
            theta = rad2deg(atan2(x(2), y(2)));
            gain_mu = FRF_data.(T(n)).grand_mean(v).gain(f);
            gain_std = FRF_data.(T(n)).grand_std(v).gain(f);
            phase_std = FRF_data.(T(n)).grand_std(v).phase(f);
            b = 2*pi*gain_mu*(phase_std/ 360);
            n_std = 1;
%             elps = Ellipse([real(mu(f)) imag(mu(f))], 2*n_std*gain_std, 2*n_std*b, ...
%                 'ellipse', theta, 0, cc(f,:), 0.4);
% 
%             elps = draw(elps);
%             set(elps.h.patch, 'EdgeColor', 'k')
%             set(elps.h.centroid, 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
%             hh.patch(v,n,f) = elps.h.patch;
%             hh.med(v,n,f) = elps.h.centroid;
        end
        if n==1
            leg(n,v) = legend(squeeze(h(v,n,:)), string(IOFv), 'Box', 'off', 'Location', 'east');
            leg(n,v).Position(1:2) = leg(n,v).Position(1:2) + [-0.2 0.04];
        end
        title(T(n))
%         uistack(squeeze(hh.patch(v,n,:)), 'top')
%         uistack(squeeze(h(v,n,:)), 'top')
    end
end

%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'Color', 'none', 'LineWidth', 1.5)

end