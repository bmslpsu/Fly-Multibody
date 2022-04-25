function [] = SOS_complex_fly()
%% SOS_complex_fly:
root = 'E:\DATA\Magno_Data\Multibody\processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select FRF data', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL')

%% FRF
clearvars -except ALL FUNC U N
clc

set_names = ["BodyFixed"];
trf_names = ["ref2head"];
yL = [1 1 1];
n_set = length(trf_names);
cond = ALL.BodyFixed.U.vel{1};
n_cond = 1;
% cc = [0.9 0 0 ; 1 0.6 0.1 ; 0.5 0.3 1 ; 0 0.4 1 ; 0 0.8 0.2 ; 0.2 0.8 1];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_cond 3*n_set])
mrksz = 6;

movegui(fig, 'center')
clear ax h hh
ax = gobjects(n_cond, n_set);
h = gobjects(n_cond, n_set, 9);
leg = gobjects(1, n_cond);
% v = 2;
for n = 1:n_set
    FRF = ALL.(set_names(n)).FRF_data.(trf_names(n));
    for v = 1
        IOFv = ALL.(set_names(n)).FRF_data.IOFv
        cmplx = FRF.fly(v).complex;
        R = real(cmplx(:));
        I = imag(cmplx(:));
        [n_freq, n_fly] = size(cmplx);
        G = repmat(IOFv, [n_fly, 1]);
        cc = hsv(n_freq);
        
        subI = v + (n-1)*n_cond;
        ax(v,n) = subplot(n_cond,n_set,subI);
        %ax(v,n) = subplot(n_plot,n_cond,subI);
        [ax_cm,h_cm] = ComplexAxes(0:0.5:yL(n));
        set(ax_cm, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
        
        h(v,n,:) = gscatter(R, I, G, cc, '.', 7, false);
        mu = FRF.grand_mean(v).complex;
        plot(real(mu), imag(mu), '-k.', 'MarkerSize', 12, 'MarkerFaceColor', 'none');
            
%         for f = 1:n_freq          
%             x = [0 real(mu(f))];
%             y = [0 real(mu(f))];
%             theta = rad2deg(atan2(x(2), y(2)));
%             gain_mu = FRF_data.(trf_names(n)).grand_mean(v).gain(f);
%             gain_std = FRF_data.(trf_names(n)).grand_std(v).gain(f);
%             phase_std = FRF_data.(trf_names(n)).grand_std(v).phase(f);
%             b = 2*pi*gain_mu*(phase_std/ 360);
%             n_std = 1;
% %             elps = Ellipse([real(mu(f)) imag(mu(f))], 2*n_std*gain_std, 2*n_std*b, ...
% %                 'ellipse', theta, 0, cc(f,:), 0.4);
% % 
% %             elps = draw(elps);
% %             set(elps.h.patch, 'EdgeColor', 'k')
% %             set(elps.h.centroid, 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
% %             hh.patch(v,n,f) = elps.h.patch;
% %             hh.med(v,n,f) = elps.h.centroid;
%         end
        if n==1
            leg(n,v) = legend(squeeze(h(v,n,:)), string(IOFv), 'Box', 'off', 'Location', 'east');
            leg(n,v).Position(1:2) = leg(n,v).Position(1:2) + [-0.2 0.04];
        end
        title(trf_names(n))
%         uistack(squeeze(hh.patch(v,n,:)), 'top')
%         uistack(squeeze(h(v,n,:)), 'top')
    end
end

%set(h.line(1:2,:),'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'Color', 'none', 'LineWidth', 1.5)

end