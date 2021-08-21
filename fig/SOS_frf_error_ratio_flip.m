function [] = SOS_frf_error_ratio_flip()
%% SOS_frf_error_ratio_flip:
root = 'E:\DATA\Magno_Data\Multibody\Processed';
[HeadFree_file,HeadFree_path] = uigetfile({'*.mat'}, ...
    'Select head free data', root, 'MultiSelect','off');

[HeadFixed_file,HeadFixed_path] = uigetfile({'*.mat'}, ...
    'Select head fixed data', root, 'MultiSelect','off');

[BodyFixed_file,BodyFixed_path] = uigetfile({'*.mat'}, ...
    'Select body fixed data', root, 'MultiSelect','off');

HeadFree = load(fullfile(HeadFree_path,HeadFree_file),'FRF_data','FUNC','U','N');
HeadFixed = load(fullfile(HeadFixed_path,HeadFixed_file),'FRF_data','FUNC','U','N');
BodyFixed = load(fullfile(BodyFixed_path,BodyFixed_file),'FRF_data','FUNC','U','N');

%% Error ratio
n_cond = HeadFree.N{1,3};
ErrorRatio = [];
for v = 1:n_cond
    % How head effects body
    ErrorRatio.fly(v).bodyfree2headfree = (HeadFree.FRF_data.ref2body.fly(v).error ./ ...
                                        HeadFree.FRF_data.ref2head.fly(v).error) - 1;
    ErrorRatio.mean(v).bodyfree2headfree = mean(ErrorRatio.fly(v).bodyfree2headfree, 2);
    ErrorRatio.std(v).bodyfree2headfree = std(ErrorRatio.fly(v).bodyfree2headfree, [], 2);
    
    % How head effects gaze
    ErrorRatio.fly(v).bodyfree2gazefree = (HeadFree.FRF_data.ref2body.fly(v).error ./ ...
                                        HeadFree.FRF_data.ref2gaze.fly(v).error) - 1;
    ErrorRatio.mean(v).bodyfree2gazefree = mean(ErrorRatio.fly(v).bodyfree2gazefree, 2);
    ErrorRatio.std(v).bodyfree2gazefree = std(ErrorRatio.fly(v).bodyfree2gazefree, [], 2);
    
    % How head-fixation effects body
    ErrorRatio.mean(v).bodyfixed2bodyfree = (HeadFree.FRF_data.ref2body.grand_mean(v).error ./ ...
                                        HeadFixed.FRF_data.ref2body.grand_mean(v).error) - 1;
    % How head-fixation effects gaze
    ErrorRatio.mean(v).bodyfixed2gazefree = (HeadFree.FRF_data.ref2gaze.grand_mean(v).error ./ ...
                                        HeadFixed.FRF_data.ref2body.grand_mean(v).error) - 1;
    % How body-fixation effects head
    ErrorRatio.mean(v).headfree2gazefixed = (HeadFree.FRF_data.ref2head.grand_mean(v).error ./ ...
                                        BodyFixed.FRF_data.ref2head.grand_mean(v).error) - 1;
    % How body-fixation effects gaze
    ErrorRatio.mean(v).gazefree2gazefixed = (HeadFree.FRF_data.ref2gaze.grand_mean(v).error ./ ...
                                        BodyFixed.FRF_data.ref2head.grand_mean(v).error) - 1;
    
end

%% All
T = string(fieldnames(ErrorRatio.mean));
n_field = length(T);
% cc = jet(n_field);
cc = distinguishable_colors(n_field);
% cc = [0 0 0; 0 0 0];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond n_field*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_field,n_cond);
pp = 1;
for n = 1:n_field
    for v = 1:n_cond
        IOFv = HeadFree.FRF_data.IOFv{v};
        ax(n,v) = subplot(n_field,n_cond,pp); hold on
            if n == 1
                title([num2str(HeadFree.U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
            end
            yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            %plot(IOFv, ErrorRatio.fly(v).(T(n)), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            [h.patch(n,v),h.line(n,v)] = PlotPatch(ErrorRatio.mean(v).(T(n)),...
                      0*ErrorRatio.mean(v).(T(n)), IOFv, 0, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
          	if v == 1
                ylabel({T(n), 'Error ratio'})
            end
            %ax(n,v).YLim(1) = -0.1;
      	pp = pp + 1;
    end
end
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
for n = 1:n_field
    linkaxes(ax(n,:), 'y')
end

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

set(ax(1,1:end),'YLim',[-1.3 1.3])
set(ax(2,1:end),'YLim',[-0.2 1])
set(ax(3:6,1:end),'YLim',1*[-1 1])

set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')

%% Head2Body error ratio
T = string(fieldnames(ErrorRatio.fly));
n_field = length(T);
% cc = jet(n_field);
% cc = distinguishable_colors(n_field, 'm');
cc = [0 0 0; 0 0 0];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond n_field*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_field,n_cond);
pp = 1;
for n = 1:n_field
    for v = 1:n_cond
        IOFv = HeadFree.FRF_data.IOFv{v};
        ax(n,v) = subplot(n_field,n_cond,pp); hold on
        if n == 1
            title([num2str(HeadFree.U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
        end
            yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            plot(IOFv, ErrorRatio.fly(v).(T(n)), '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            [h.patch(n,v),h.line(n,v)] = PlotPatch(ErrorRatio.mean(v).(T(n)),...
                      ErrorRatio.std(v).(T(n)), IOFv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
          	if v == 1
                ylabel({T(n), 'Error ratio'})
            end
            
      	pp = pp + 1;
    end
end
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
for n = 1:n_field
    linkaxes(ax(n,:))
end

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

set(ax(1,1:end),'YLim',[-1.5 1.5])
set(ax(2,1:end),'YLim',[-0.2 1])

set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')

%% Head-fixed body to head-free body & gaze
T = ["bodyfixed2bodyfree", "bodyfixed2gazefree"];
n_field = length(T);
% cc = distinguishable_colors(n_field, 'm');
% cc = [0 0.7 0.7; 0.9 0.9 0];
cc = [0 0 0; 0.4 0.4 0.6];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond 1*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_cond);
for n = 1:n_field
    for v = 1:n_cond
        IOFv = HeadFree.FRF_data.IOFv{v};
        ax(1,v) = subplot(1,n_cond,v); hold on
            %title([num2str(HeadFree.U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
            yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            %plot(IOFv, ErrorRatio.fly(v).(T(n)), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            [h.patch(n,v),h.line(n,v)] = PlotPatch(ErrorRatio.mean(v).(T(n)),...
                      0*ErrorRatio.mean(v).(T(n)), IOFv, 0, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
          	if v == 1
                ylabel({T(n), 'Error ratio'})
            end  
    end
    ax(1,v).YLim(1) = -0.1;
end
leg = legend(h.line(:,1), T, 'Box', 'off');
leg.Position = [0.38 0.73 0.23 0.18];
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'xy')

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

set(ax(1,1:end),'YLim',[-1 1])

set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')

%% Body-fixed head to head-free head & gaze
T = ["headfree2gazefixed", "gazefree2gazefixed"];
n_field = length(T);
% cc = prism(n_field);
% cc = [0 0.7 0.7; 0.9 0.9 0];
cc = cool(n_field);
cc = [0 0 0; 0.4 0.4 0.6];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_cond 1*1.8])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_cond);
for n = 1:n_field
    for v = 1:n_cond
        IOFv = HeadFree.FRF_data.IOFv{v};
        ax(1,v) = subplot(1,n_cond,v); hold on
            %title([num2str(HeadFree.U{1,3}{1}(v)) '°/s'], 'interpreter', 'none')
            yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            %plot(IOFv, ErrorRatio.fly(v).(T(n)), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            [h.patch(n,v),h.line(n,v)] = PlotPatch(ErrorRatio.mean(v).(T(n)),...
                      0*ErrorRatio.mean(v).(T(n)), IOFv, 0, 1, cc(n,:), 0.7*cc(n,:), 0.2, 1);
                  
          	if v == 1
                ylabel('Error ratio')
            end  
    end
    ax(1,v).YLim(1) = -0.1;
end
leg = legend(h.line(:,1), T, 'Box', 'off');
leg.Position = [0.66 0.76 0.23 0.18];
set(h.line, 'Marker', '.','MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 0.75)
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'XLim', [0.2 20],...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')
set(ax, 'XTick', [0.1, 1 10])

linkaxes(ax, 'x')
for n = 1:n_field
    linkaxes(ax(1,:))
end

XLabelHC = get(ax(end,:), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')

set(ax(1,1:end),'YLim',[-1 1])

set(ax(:,2:end), 'YColor', 'none')

set(ax,'XScale','log')

% Head fixation causes an increase in performance at low frequencies at the detriment/conseqeunce of
% a decrease in performance at high frequencies

%% Save time constant data
% head_free_filedata = textscan(HeadFree_file, '%s', 'delimiter', '._');
% head_fixed_filedata = textscan(HeadFixed_file, '%s', 'delimiter', '._');
% body_fixed_filedata = textscan(BodyFixed_file, '%s', 'delimiter', '._');
% head_free_name = [];
% head_fixed_name = [];
% body_fixed_name = [];
% for n = 2:5
%     head_free_name = [head_free_name '_' char(head_free_filedata{1}(n))];
%     head_fixed_name = [head_fixed_name '_' char(head_fixed_filedata{1}(n))]; 
%     body_fixed_name = [body_fixed_name '_' char(body_fixed_filedata{1}(n))]; 
% end
% fname = ['TimeConstant' head_free_name head_fixed_name body_fixed_name]; 
% 
% root = 'E:\DATA\Magno_Data\Multibody';
% savedir = fullfile(root,'processed');
% mkdir(savedir)
% save(fullfile(savedir, [fname '.mat']), 'time_const', 'time_const_all', 'r2', 'r2_all');
end