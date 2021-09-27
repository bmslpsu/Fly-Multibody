function [] = SOS_head_free2fixed()
%% SOS_head_free2fixed:
%
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'ALL');

%% Make stats table
clearvars -except PATH FILE ALL
clc

cond = ALL.HeadFree.U{1,3}{1};
n_cond = ALL.HeadFree.N{1,3};


clss = ["BodyFixed", "HeadFree"];
trf = ["ref2head", "ref2head"];

n_cond = ALL.HeadFree.N{1,3};
IOFv = ALL.(clss(1)).FRF_data.IOFv;

T = [];
for v = 1:n_cond
    T.grand_mean(v).complex = ALL.(clss(1)).FRF_data.(trf(1)).grand_mean(v).complex ./ ...
        ALL.(clss(2)).FRF_data.(trf(2)).grand_mean(v).complex ;
    T.grand_std(v).complex = sqrt( ALL.(clss(1)).FRF_data.(trf(1)).grand_std(v).complex.^(2) + ...
        ALL.(clss(2)).FRF_data.(trf(2)).grand_std(v).complex.^(2) ) ;
    
    T.grand_mean(v).gain = abs(T.grand_mean(v).complex);
    T.grand_std(v).gain = abs(T.grand_std(v).complex);
    
    T.grand_mean(v).phase = rad2deg(angle(T.grand_mean(v).complex));
    T.grand_std(v).phase = rad2deg(angle(T.grand_std(v).complex));
end
% free2fixed.ref2head = T;


plot_names = ["gain", "phase"];
n_plot = length(plot_names);
cc = [0.6 0.1 0.4];

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches')
movegui(fig, 'center')
clear ax h
ax = gobjects(n_plot,n_cond);
for v = 1:n_cond
    subI = v + (0:n_plot-1)*n_cond;
    for p = 1:n_plot
        ax(p,v) = subplot(n_plot,n_cond,subI(p)); cla ; hold on
        %title([num2str(cond(v)) '°/s'], 'interpreter', 'none')
            [h.patch(p,v),h.line(p,v)] = PlotPatch(...
                    T.grand_mean(v).(plot_names(p)),...
                  	T.grand_std(v).(plot_names(p)), ...
                    IOFv{v}, 1, 1, cc, 0.7*cc, 0.2, 1);
    end
end
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1)
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 10, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])
set(ax,'XScale','log')

set(ax(1:end-1,:), 'XColor', 'none')
set(ax(:,2:end), 'YColor', 'none')

end