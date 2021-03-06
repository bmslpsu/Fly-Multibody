function [] = make_head_free_magno_SOS_vel(rootdir)
%% MakeData_SOS_v1_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%
rootdir = 'E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v1';
filename = 'SOS_HeadFree';

%% Setup Directories %%
root.daq = rootdir; clear rootdir
root.body = fullfile(root.daq,'tracked_body');
root.reg = fullfile(root.daq,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head_tip');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root.head,'*.mat',false);
% [D,I,N,U,T,~,~,basename] = GetFileData(root.benifly,'*.csv',false);

% Load function files
FUNC = cell(N.vel, 1);
FUNC{1} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v1\function\ALL_position_function_SOS_Fs_50_T_20_vel_62_freq_8.2_5.05_3.1_1.9_1.15_0.7_0.45_0.25_0.15_amp_1.2_1.95_3.2_5.25_8.55_14_22.85_37.35_61.mat');
FUNC{2} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v1\function\ALL_position_function_SOS_Fs_50_T_20_vel_103_freq_9.1_6.1_4.05_2.7_1.8_1.2_0.8_0.55_0.35_amp_1.8_2.7_4_6_9_13.45_20.1_30.1_45.mat');
FUNC{3} = load('E:\EXPERIMENTS\MAGNO\Experiment_SOS_vel_v1\function\ALL_position_function_SOS_Fs_50_T_20_vel_148_freq_9.4_6.3_4.25_2.85_1.9_1.3_0.85_0.6_0.4_amp_2.5_3.75_5.55_8.3_12.35_18.4_27.45_40.9_61.mat');

%% Get Data %%
close all
clc

Fs = 100;
Fc = 20;
func_length = 20;
tintrp = (0:(1/Fs):func_length)';
debug = false;
[b,a] = butter(3, Fc/(Fs/2),'low');
ALL = cell(N.fly,N.vel);
DATA = [I , splitvars(table(num2cell(zeros(N.file,5))))]; % store saccade objects
DATA.Properties.VariableNames(4:end) = {'reference','body','head','error','dwba',};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    % Load DAQ, body, head, & wing data
	data.daq = load(fullfile(root.daq,  [basename{n} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body = load(fullfile(root.body, [basename{n} '.mat']),'bAngles'); % load body angles
	data.head = load(fullfile(root.head, [basename{n} '.mat']),'head_data'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...
                            [basename{n} '.csv'])); % load head & wing angles from Benifly
    
    % Get synced frame times and pattern data
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, func_length,trigger, true, [], false, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    lwing = hampel(data.benifly.Time, data.benifly.LWing);
    rwing = hampel(data.benifly.Time, data.benifly.RWing);
	lwing = rad2deg(filtfilt(b,a,lwing));
    rwing = rad2deg(filtfilt(b,a,rwing));
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.bAngles;
	head = data.head.head_data.angle;
    
    % Interpolate so all signals have the same times
    [b_pat, a_pat] = butter(3, 12 / (Fs/2), 'low');
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = filtfilt(b_pat, a_pat, Reference);
    Reference = Reference - mean(Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = filtfilt(b, a, Head);
    Error   = Reference - Body - Head;
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    % Store signals
    n_detrend = 1;
    DATA.reference{n}  = singal_attributes(Reference, tintrp);
    DATA.body{n}       = singal_attributes(Body, tintrp, 20, n_detrend);
    DATA.head{n}       = singal_attributes(Head, tintrp, 20);
    DATA.error{n}      = singal_attributes(Error, tintrp, 20, n_detrend);
    DATA.dwba{n}       = singal_attributes(dWBA, tintrp, 20);

    % Debug plot.
    if debug
        figure (100)
        clear ax
        ax(1) = subplot(1,1,1) ; cla ; hold on
            plot(tintrp, DATA.reference{n}.position, 'k', 'LineWidth', 1)
            %plot(tintrp, DATA.error{kk}.detrend, 'g', 'LineWidth', 1)
            plot(tintrp, DATA.body{n}.position, 'r', 'LineWidth', 1)
            plot(tintrp, DATA.head{n}.position, 'b', 'LineWidth', 1)
            plot(tintrp, 5*DATA.dwba{n}.position_lpf, 'm', 'LineWidth', 1)
            leg = legend('Reference','Body','Head','\DeltaWBA', 'Orientation', 'horizontal');
            xlabel('Time (s)')
            ylabel('(°)')
                   
       	set(gcf, 'Color', 'w')
       	set(ax, 'Linewidth', 2)
        linkaxes(ax,'x')
        pause
    end
    
    IOFreq = flipud(FUNC{I.vel(n)}.All.Freq);
    SYS_ref2_head_body  = frf(tintrp, Reference, IOFreq, false, Body, Head);
    %pause
    %close all
    SYS_ref2_wing       = frf(tintrp, Reference, IOFreq, false, LWing, RWing, dWBA);
    SYS_head2_body_wing = frf(tintrp, Head, IOFreq, false, Body, dWBA);
    SYS_wing2_body      = frf(tintrp, dWBA, IOFreq, false, Body);
    
    SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_ref2_wing, ...
                                    SYS_head2_body_wing, SYS_wing2_body);
    
    ALL{I.fly(n),I.vel(n)}(end+1,1) = SYS_all;
end

%% Group Data
clc
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
for v = 1:N.vel
    GRAND.all(v) = cell2struct(cell(nfield,1),fields);
    GRAND.all_trial(v) = cell2struct(cell(nfield,1),fields);
    %GRAND.fly_all(fr) = cell2struct(cell(n_stat_fields,1),stat_fields);
    %GRAND.fly_all(fr) = [];
    for n = 1:N.fly
        for f = 1:nfield
            FLY.all(n,v).(fields{f})   = cat(3,ALL{n,v}.(fields{f}));
            FLY.stats(n,v).(fields{f}) = system_stats(FLY.all(n,v).(fields{f}),3);
            GRAND.all(v).(fields{f}) 	= cat(3,GRAND.all(v).(fields{f}), ...
                                                        FLY.all(n,v).(fields{f}));

            stat_fields = fieldnames(FLY.stats(n,v).(fields{f}));
            n_stat_fields = length(stat_fields);
            for s = 1:n_stat_fields
                GRAND.fly_all(v).(stat_fields{s}).(fields{f})(:,:,n) = ...
                    FLY.stats(n,v).(fields{f}).(stat_fields{s});
                
                GRAND.fly_stats(v).(stat_fields{s}).(fields{f}) = ...
                    system_stats(GRAND.fly_all(v).(stat_fields{s}).(fields{f}),3);
            end
        end
    end
    GRAND.all_trial(v) = structfun(@(x) system_stats(x,3), GRAND.all(v), 'UniformOutput', false);
end

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Magno_Data\Multibody';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'FUNC', 'DATA', 'ALL', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')

%%
clc
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 9 4])
movegui(fig, 'center')
clear ax h
axis tight
xlabel('Time (s)')
ylabel('Pattern position (°)')
% x_start = All.X_step(1);
%sync = DATA.reference{n}.position(1);
% plot(All.time, All.X_step, 'k', 'LineWidth', 1)

ax = gobjects(N.vel,1);
for v = 1:N.vel
    ax(v) = subplot(N.vel,1,v) ; hold on
    tt = FUNC{v}.All.time;
    pos = FUNC{v}.All.X;
    plot(tt, pos, 'k', 'LineWidth', 1)
end

test = cell(N.vel,1);
for n = 1:N.file
    vel = DATA.vel(n);
    DATA.reference{n}.position(1) = DATA.reference{n}.position(2);
    
    ax(vel) = subplot(N.vel,1,vel) ; hold on
    plot(DATA.reference{n}.time, DATA.reference{n}.position, 'LineWidth', 1);
    %test{vel}(:,end+1) = DATA.reference{n}.position;
    test{vel}(:,end+1) = DATA.body{n}.mag.position;
    %test{vel}(:,end+1) = DATA.dwba{n}.mag.velocity;
end


% legend('Computer generated', 'Experimental (from panel controller)', 'Box', 'off')
Fv = DATA.reference{1}.Fv;
% plot(Fv, test)
set(ax, 'LineWidth', 1.5, 'XLim', [-0.2 20.2])

%%
bI = 1;
hI = 2;
gI = 3;
lwI = 4;
rwI = 5;
dwI = 6;
h2b = 8;
h2w = 9;
w2b = 11;

rI = 1;
hcolor = 'b';
bcolor = 'r';
gcolor = [0.4 0.1 0.7];
ccolor = 'k';
rcolor = 'g';

%% FRF
pI = [bI hI gI dwI h2b h2w w2b];
T = ["ref2body", "ref2head", "ref2gaze", "ref2wing", "head2body", "head2wing", "wing2body"];
n_plot = length(pI);
cc = hsv(n_plot);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*n_plot 8])
movegui(fig, 'center')
clear ax h
ax = gobjects(3,n_plot);
phase_lim = [nan nan nan 50 nan 50 0];
shift_I = {7, 7, 7, nan, 7, nan, 7};
for n = 1:n_plot
    subI = n + (0:2)*n_plot;
    ax(1,n) = subplot(3,n_plot,subI(1)); hold on ; title(T(n))
        [h.patch(1,n),h.line(1,n)] = PlotPatch(grand_mean_gain(:,pI(n)),...
                  grand_std_gain(:,pI(n)), Fv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 2);

    ax(2,n) = subplot(3,n_plot,subI(2)); hold on
        plot([0 10], [0 0], '--k')
        phase = grand_mean_phase(:,pI(n));
        phase(phase > phase_lim(n)) = phase(phase > phase_lim(n)) - 360;
        phase(any((1:N.freq)'==shift_I{n},2)) = phase(any((1:N.freq)'==shift_I{n},2)) - 360;
        [h.patch(2,n),h.line(2,n)] = PlotPatch(phase,...
                  grand_std_phase(:,pI(n)), Fv, 1, 1,cc(n,:), 0.7*cc(n,:), 0.2, 2);

    ax(3,n) = subplot(3,n_plot,subI(3)); hold on
        [h.patch(3,n),h.line(3,n)] = PlotPatch(grand_mean_cohr(:,pI(n)),...
                  grand_std_cohr(:,pI(n)), Fv, 1, 1, cc(n,:), 0.7*cc(n,:), 0.2, 2);

end
set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontSize', 12, 'XLim', [0 5.4],...
    'XGrid','on','YGrid','on','box','on')

linkaxes(ax,'x')
% linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
linkaxes(ax(3,:),'y')

YLabelHC = get(ax(3,:), 'XLabel');
set([YLabelHC{:}], 'String', 'Frequency (Hz)')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Gain (°/°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Phase Difference (°)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Coherence')

% set(ax(1:3),'YLim',[0 1.05])
% set(ax(2), 'YLim', [0 0.5])
% set(ax(4:6),'YLim',[-360 120])
% set(ax(7:9),'YLim',[0 1.05])
% set(ax(1:4),'XTickLabels',[])

set(ax,'XScale','log')
align_Ylabels(fig)

%% Complex Gain
clc
grand_stats_all = cat(1, GRAND.all_trial.IOFRF);
grand_stats_all = cat(1, grand_stats_all.mean);
pI = [bI hI gI dwI h2b h2w w2b];
T = ["ref2body", "ref2head", "ref2gaze", "ref2wing", "head2body", "head2wing", "wing2body"];
n_plot = length(pI);

clear ax h
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3*ceil(n_plot/2) 6])
ax = gobjects(1,n_plot);
h = gobjects(N.freq,n_plot);
cc = hsv(N.freq);
for n = 1:n_plot
    for v = 1:N.freq
        all_comp = cat(3, GRAND.all(v).IOFRF);
        ax(1,n) = subplot(2,ceil(n_plot/2),n);
            all = squeeze(all_comp(:,pI(n),:));
            %[~,~] = ComplexAxes(0:0.2:1.2);
%             h(fr,n) = plot(real(all), imag(all), '.', 'MarkerSize', 10, ...
%                 'MarkerFaceColor', 'none', 'MarkerEdgeColor', cc(fr,:));
           h(v,n) = polarplot(angle(all), abs(all), '.', 'MarkerSize', 10, ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', 0.3*cc(v,:));
            hold on ; title(T(n))
            polarplot(angle(grand_stats_all(v,pI(n))), abs(grand_stats_all(v,pI(n))), ...
                '.', 'MarkerSize', 25, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', cc(v,:));
            hold on ; title(T(n))
            
%         h.mean(1) = plot(real(grand_stats_all(:,bI)), imag(grand_stats_all(:,bI)), ...
%                     '.', 'MarkerSize', 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
    end
end
% leg = legend(h(:,end), string(U.freq{1}), 'Box', 'on', 'Location', 'east');
% leg.Title.String = 'Frequency (Hz)';
set(ax, 'LineWidth', 1.5)

%% Time Domain
time = cat(1, GRAND.all_trial.Time);
time = cat(2, time.mean);
time = mean(time,2);
grand_stats_all = cat(1, GRAND.all_trial.State);
grand_stats_mean = cat(3, grand_stats_all.mean);
grand_stats_std = cat(3, grand_stats_all.std);
head = squeeze(grand_stats_mean(:,hI,:));
body = squeeze(grand_stats_mean(:,bI,:));
head_std = squeeze(grand_stats_std(:,hI,:));
body_std = squeeze(grand_stats_std(:,bI,:));
ref = cat(1, GRAND.all_trial.refState);
ref = cat(3, ref.mean);
ref = squeeze(ref(:,1,:));
Amp = flipud(U.Amp{1});

clear ax h
cc = hsv(N.freq);
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 10 9])
movegui(fig, 'center')
ax = gobjects(N.freq,1);
for v = 1:N.freq
    T = ['Freq = ' num2str(U.freq{1}(v)) 'Hz   Amp = ' num2str(Amp(v)) '°'];
    ax(v) = subplot(N.freq,1,v); hold on ; title(T)
        plot(time, ref(:,v), 'k', 'LineWidth', 1)
        [h.patch(v),h.line(v)] = PlotPatch(body(:,v),...
          body_std(:,v), time, 1, 1, bcolor, bcolor, 0.2, 1);
end
set(ax, 'LineWidth', 1.5)
linkaxes(ax, 'x')

%% Weights
clear ax h
fig(4) = figure (4) ; clf
set(fig(4),'Color','w','Units','inches','Position',[4 2 4 4])
ax(1) = subplot(1,1,1); hold on ; xlabel('Frequency (Hz)') ; ylabel('Weights')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOWeight(:,bI,:)), ...
        'Color', [1 0 0 0.3])
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOWeight(:,hI,:)), ...
        'Color', [0 0 1 0.3])

  	[h.patch(1), h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOWeight.mean(:,bI),...
              GRAND.fly_stats.std.IOWeight.mean(:,bI), GRAND.fly_stats.mean.IOFv.mean(:,1), ...
              1, 1, bcolor,'r', 0.2, 2);
  	[h.patch(2), h.line(2)] = PlotPatch(GRAND.fly_stats.mean.IOWeight.mean(:,hI),...
              GRAND.fly_stats.std.IOWeight.mean(:,hI), GRAND.fly_stats.mean.IOFv.mean(:,1), ...
              1, 1, hcolor,'b', 0.2, 2);
          
set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax, 'YLim', [0 1])
legend(h.line, 'Body', 'Head', 'Box', 'off')

%% Time
clear ax h
fig(5) = figure (5) ; clf
set(fig(5),'Color','w','Units','inches','Position',[2 2 8 5])
ax(1) = subplot(3,1,1); hold on ; ylabel('Body (°)')
    plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.State(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
 	plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.refState(:,rI,:)), ...
        'Color', rcolor, 'LineWidth', 1)
  	[h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.State.mean(:,bI),...
              GRAND.fly_stats.std.State.mean(:,bI), GRAND.fly_stats.mean.Time.mean(:,1), ...
              1, 1, bcolor,'b', 0.2, 2);
	uistack(h.patch(1), 'bottom')

ax(2) = subplot(3,1,2); hold on ; ylabel('Head (°)') ; ylim(15.5*[-1 1]) ; yticks(15*[-1 1])
    plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.State(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
 	plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.refState(:,rI,:)), ...
        'Color', rcolor, 'LineWidth', 1)
  	[h.patch(2),h.line(2)] = PlotPatch(GRAND.fly_stats.mean.State.mean(:,hI),...
              GRAND.fly_stats.std.State.mean(:,hI), GRAND.fly_stats.mean.Time.mean(:,1), ...
              1, 1, hcolor,'b', 0.2, 2);
	uistack(h.patch(2), 'bottom')
          
ax(3) = subplot(3,1,3); hold on ; ylabel('Gaze (°)') ; xlabel('Time (s)')
    plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.State(:,gI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
 	plot(squeeze(GRAND.all.Time(:,1,:)), squeeze(GRAND.all.refState(:,rI,:)), ...
        'Color', rcolor, 'LineWidth', 1)
  	[h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.State.mean(:,gI),...
              GRAND.fly_stats.std.State.mean(:,gI), GRAND.fly_stats.mean.Time.mean(:,1), ...
              1, 1, gcolor,'b', 0.2, 2);
    uistack(h.patch(3), 'bottom')
    
set(ax,'LineWidth',1.5,'FontWeight','Bold','FontSize',10,'XLim',[0 20])
linkaxes(ax,'x')
linkaxes(ax([1,3]),'y')
set(ax(1), 'YLim', 75*[-1 1])

%% Wing Complex Gain
clear ax
cc = hsv(length(IOFreq));
fig(6) = figure (2) ; clf
set(fig(6), 'Color', 'w', 'Units', 'inches', 'Position', [2 2 10 5])
ax(1) = subplot(1,3,1); hold on ; title('Left Wing')
    Real = real(squeeze(GRAND.all.IOFRF(:,lwI,:)))';
    Img  = imag(squeeze(GRAND.all.IOFRF(:,lwI,:)))';
    [~,~] = ComplexAxes(0:0.05:0.2);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
ax(2) = subplot(1,3,2); hold on ; title('Right Wing')
    Real = real(squeeze(GRAND.all.IOFRF(:,rwI,:))');
    Img  = -imag(squeeze(GRAND.all.IOFRF(:,rwI,:))');
    [~,~] = ComplexAxes(0:0.05:0.2);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
ax(3) = subplot(1,3,3); hold on ; title(' \DeltaWBA')
    Real = real(squeeze(GRAND.all.IOFRF(:,rwI,:))');
    Img  = -imag(squeeze(GRAND.all.IOFRF(:,dwI,:))');
    [~,~] = ComplexAxes(0:0.05:0.2);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
    
set(ax, 'LineWidth',1.5, 'FontWeight', 'Bold', 'FontSize', 12)

%% Wing FRF
clear ax h
fig(7) = figure (7) ; clf
set(fig(7), 'Color', 'w', 'Units', 'inches', 'Position', [4 1 3*3 3*2])
ax(1) = subplot(3,3,1); hold on ; ylabel('Gain (°/°)') ; title('Left Wing')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,lwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,lwI),...
              GRAND.fly_stats.std.IOGain.mean(:,lwI), IOFreq, ...
              1, 1, bcolor, 'b', 0.2, 2);
          
ax(2) = subplot(3,3,2); hold on ; title('Right Wing')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,rwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,rwI),...
              GRAND.fly_stats.std.IOGain.mean(:,rwI), IOFreq, ...
              1, 1, hcolor, 'b', 0.2, 2);
          
ax(3) = subplot(3,3,3); hold on ; title('\DeltaWBA')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,dwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,dwI),...
              GRAND.fly_stats.std.IOGain.mean(:,dwI), IOFreq, ...
              1, 1, gcolor, 'b', 0.2, 2);       
                 
ax(4) = subplot(3,3,4); hold on ; ylabel('Phase Difference (°)')
    plot([0 10], [0 0], '--k')
    phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,lwI));
    phase(phase > 100) = phase(phase > 100) - 360;
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,lwI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(4),h.line(4)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,lwI)), IOFreq, ...
              1, 1, bcolor, 'b', 0.2, 2);
          
ax(5) = subplot(3,3,5); hold on
    plot([0 10], [0 0], '--k')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,rwI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    
	phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,rwI));
    %(phase>0) = phase(phase>0) - 360;
    [h.patch(5),h.line(5)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,rwI)), IOFreq, ...
              1, 1, hcolor, 'b', 0.2, 2);
          
ax(6) = subplot(3,3,6); hold on
    plot([0 10], [0 0], '--k')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,dwI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    
	phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,dwI));
    %(phase>0) = phase(phase>0) - 360;
    [h.patch(6),h.line(6)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,dwI)), IOFreq, ...
              1, 1, gcolor, 'b', 0.2, 2);
               
ax(7) = subplot(3,3,7); hold on ; ylabel('Coherence')
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,lwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,lwI),...
              GRAND.fly_stats.std.Cohr.mean(:,lwI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, bcolor,'b', 0.2, 1);
    [h.patch(7),h.line(7)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,lwI),...
              GRAND.fly_stats.std.IOCohr.mean(:,lwI), IOFreq, ...
              1, 1, ccolor, 'b', 0.2, 2);
          
ax(8) = subplot(3,3,8); hold on
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,rwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,rwI),...
              GRAND.fly_stats.std.Cohr.mean(:,rwI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, hcolor,'b', 0.2, 1);
    [h.patch(8),h.line(8)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,rwI),...
              GRAND.fly_stats.std.IOCohr.mean(:,rwI), IOFreq, ...
              1, 1, ccolor, 'b', 0.2, 2);

ax(9) = subplot(3,3,9); hold on
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,dwI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,dwI),...
              GRAND.fly_stats.std.Cohr.mean(:,dwI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, gcolor, 'b', 0.2, 1);
    [h.patch(9),h.line(9)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,dwI),...
              GRAND.fly_stats.std.IOCohr.mean(:,dwI), IOFreq, ...
              1, 1, ccolor, 'b', 0.2, 2);
          
set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax(1:3),'YLim',[0 0.3])
set(ax(4:6),'YLim',[-300 150])
set(ax(7:9),'YLim',[0 1])
linkaxes(ax,'x')
linkaxes(ax(1:3),'y')
linkaxes(ax(4:6),'y')
linkaxes(ax(7:9),'y')
% set(ax(1:4),'XTickLabels',[])
YLabelHC = get(ax(7:9), 'XLabel');
set([YLabelHC{:}], 'String', 'Frequency (Hz)')
% set(ax,'XScale','log')
align_Ylabels(fig(7))

%% Head 2 Wing FRF
h2w_color = 'c';
clear ax h
fig(8) = figure (8) ; clf
set(fig(8), 'Color', 'w', 'Units', 'inches', 'Position', [4 1 3 3*2])
ax(1) = subplot(3,1,1); hold on ; ylabel('Gain (°/°)') ; title('Head >>> Wing')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,h2w,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,h2w),...
              GRAND.fly_stats.std.IOGain.mean(:,h2w), IOFreq, ...
              1, 1, h2w_color, 'b', 0.2, 2);
                        
ax(2) = subplot(3,1,2); hold on ; ylabel('Phase Difference (°)')
    plot([0 10], [0 0], '--k')
    phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,h2w));
    phase(phase > 100) = phase(phase > 100) - 360;
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,h2w,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,h2w)), IOFreq, ...
              1, 1, h2w_color, 'b', 0.2, 2);
          
ax(3) = subplot(3,1,3); hold on ; ylabel('Coherence')
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,h2w,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,h2w),...
              GRAND.fly_stats.std.Cohr.mean(:,h2w), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, h2w_color, 'b', 0.2, 1);
%     [h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,h2w),...
%               GRAND.fly_stats.std.IOCohr.mean(:,h2w), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);

set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')  
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax(2),'YLim',[-200 100])
align_Ylabels(fig(8))

%% Head 2 Body FRF
h2b_color = 'm';
clear ax h
fig(9) = figure (9) ; clf
set(fig(9), 'Color', 'w', 'Units', 'inches', 'Position', [4 1 3 3*2])
ax(1) = subplot(3,1,1); hold on ; ylabel('Gain (°/°)') ; title('Head >>> Body')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,h2b,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,h2b),...
              GRAND.fly_stats.std.IOGain.mean(:,h2b), IOFreq, ...
              1, 1, h2b_color, 'b', 0.2, 2);
                        
ax(2) = subplot(3,1,2); hold on ; ylabel('Phase Difference (°)')
    plot([0 10], [0 0], '--k')
    phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,h2b));
    phase(phase > 100) = phase(phase > 100) - 360;
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,h2b,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,h2b)), IOFreq, ...
              1, 1, h2b_color, 'b', 0.2, 2);
          
ax(3) = subplot(3,1,3); hold on ; ylabel('Coherence')
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,h2b,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,h2b),...
              GRAND.fly_stats.std.Cohr.mean(:,h2b), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, h2b_color, 'b', 0.2, 1);
%     [h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,h2w),...
%               GRAND.fly_stats.std.IOCohr.mean(:,h2w), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);

set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax(2),'YLim',[-200 100])
align_Ylabels(fig(9))

%% Wing 2 Body FRF
w2b_color = [0.2 0.8 0.4];
clear ax h
fig(10) = figure (10) ; clf
set(fig(10), 'Color', 'w', 'Units', 'inches', 'Position', [4 1 3 3*2])
ax(1) = subplot(3,1,1); hold on ; ylabel('Gain (°/°)') ; title('Wing >>> Body')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,w2b,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,w2b),...
              GRAND.fly_stats.std.IOGain.mean(:,w2b), IOFreq, ...
              1, 1, w2b_color, 'b', 0.2, 2);
                        
ax(2) = subplot(3,1,2); hold on ; ylabel('Phase Difference (°)')
    plot([0 10], [0 0], '--k')
    phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,w2b));
    phase(phase > 100) = phase(phase > 100) - 360;
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,w2b,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,w2b)), IOFreq, ...
              1, 1, w2b_color, 'b', 0.2, 2);
          
ax(3) = subplot(3,1,3); hold on ; ylabel('Coherence')
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,w2b,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,w2b),...
              GRAND.fly_stats.std.Cohr.mean(:,w2b), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, w2b_color, 'b', 0.2, 1);
%     [h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,h2w),...
%               GRAND.fly_stats.std.IOCohr.mean(:,h2w), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);

set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax(2),'YLim',[-200 100])
align_Ylabels(fig(10))

%% SAVE %%
% disp('Saving...')
% save(['H:\DATA\Magno_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'], ...
%     'FLY','GRAND','D','I','U','N','T','-v7.3')
% disp('SAVING DONE')
end