function [] = MakeData_SOS_v1_HeadFree(rootdir)
%% MakeData_SOS_v1_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%
rootdir = 'E:\Experiment_SOS_v1';
filename = 'SOS_HeadFree_DATA';

%% Setup Directories %%
root.daq = rootdir; clear rootdir
root.body = fullfile(root.daq,'tracked_body');
root.reg = fullfile(root.daq,'registered');
root.benifly = fullfile(root.reg ,'tracked_head_wing');
root.head = fullfile(root.reg ,'tracked_head');

% Select files
[FILES, ~] = uigetfile({'*.csv'}, 'Select trials', root.benifly, 'MultiSelect','on');
FILES = cellstr(FILES)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,'',false);

%% Get Data %%
close all
clc

IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
Fs = 100;
Fc = 15;
tintrp = (0:(1/Fs):20)';  
function_length = 20;
[b,a] = butter(2, Fc/(Fs/2),'low');
ALL = cell(N.fly,1);
debug = false;
for kk = 1:N{1,end}
    disp(kk)
    % Load DAQ, body, head, & wing data
	data.daq     = load(fullfile(root.daq,  [basename{kk} '.mat']),'data','t_p'); % load camera trigger & pattern x-position
    data.body    = load(fullfile(root.body, [basename{kk} '.mat']),'bAngles'); % load body angles
	data.head    = load(fullfile(root.head, [basename{kk} '.mat']),'hAngles'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...
                            [basename{kk} '.csv'])); % load head & wing angles from Benifly
    
    % Sync all the signals
    daq_time    = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger     = data.daq.data(:,1);
    
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, function_length, trigger, tintrp, [], false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    lwing = hampel(data.benifly.Time, data.benifly.LWing);
    rwing = hampel(data.benifly.Time, data.benifly.RWing);
	lwing = rad2deg(filtfilt(b,a,lwing));
    rwing = rad2deg(filtfilt(b,a,rwing));
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    body = data.body.bAngles;
	head = data.head.hAngles;
    % head = rad2deg(-data.benifly.Head);
    
    % Interpolate so all signals have the same times
    Reference = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Reference = Reference - mean(Reference);
    
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Body    = Body - mean(Body);
    
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    Head    = Head - mean(Head);
    
    Error   = Reference - Body - Head;
    
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    
    dWBA    = -interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    dWBA    = dWBA - mean(dWBA);
    
    SYS_ref2_head_body  = frf(tintrp, Reference, IOFreq, debug, Body, Head);
    SYS_ref2_wing       = frf(tintrp, Reference, IOFreq, debug, LWing, RWing, dWBA);
    SYS_head2_body_wing = frf(tintrp, Head, IOFreq, debug, Body, dWBA);
    SYS_wing2_body      = frf(tintrp, dWBA, IOFreq, debug, Body);
    
    SYS_all = CatStructFields(2, SYS_ref2_head_body, SYS_ref2_wing, SYS_head2_body_wing,SYS_wing2_body);
    
    ALL{I.fly(kk)}(end+1,1) = SYS_all;
end

%% Group Data
fields = fieldnames(ALL{1});
nfield = length(fields);
FLY = [];
GRAND = [];
GRAND.all = cell2struct(cell(nfield,1),fields);
for kk = 1:N.fly
    for f = 1:nfield
        FLY.all(kk).(fields{f})   = cat(3,ALL{kk}.(fields{f}));
        FLY.stats(kk).(fields{f}) = system_stats(FLY.all(kk).(fields{f}),3);
        GRAND.all.(fields{f})     = cat(3,GRAND.all.(fields{f}), FLY.all(kk).(fields{f}));
       
        stat_fields = fieldnames(FLY.stats(kk).(fields{f}));
        n_stat_feilds = length(stat_fields);
        for s = 1:n_stat_feilds
            GRAND.fly_all.(stat_fields{s}).(fields{f})(:,:,kk) = FLY.stats(kk).(fields{f}).(stat_fields{s});
            GRAND.fly_stats.(stat_fields{s}).(fields{f}) = ...
                system_stats(GRAND.fly_all.(stat_fields{s}).(fields{f}),3);
        end
    end
end
GRAND.all_trial = structfun(@(x) system_stats(x,3), GRAND.all, 'UniformOutput', false);

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

%% Complex Gain
clear ax
cc = hsv(length(IOFreq));
fig(2) = figure (2) ; clf
set(fig(2),'Color','w','Units','inches','Position',[2 2 10 5])
ax(1) = subplot(1,3,1); hold on ; title('Body')
    Real = real(squeeze(GRAND.all.IOFRF(:,bI,:)))';
    Img  = imag(squeeze(GRAND.all.IOFRF(:,bI,:)))';
    [~,~] = ComplexAxes(0:0.2:1.2);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
ax(2) = subplot(1,3,2); hold on ; title('Head')
    Real = real(squeeze(GRAND.all.IOFRF(:,hI,:))');
    Img  = -imag(squeeze(GRAND.all.IOFRF(:,hI,:))');
    [~,~] = ComplexAxes(0:0.2:0.5);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
ax(3) = subplot(1,3,3); hold on ; title('Gaze')
    Real = real(squeeze(GRAND.all.IOFRF(:,gI,:))');
    Img  = -imag(squeeze(GRAND.all.IOFRF(:,gI,:))');
    [~,~] = ComplexAxes(0:0.2:1);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    %legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
set(ax, 'LineWidth',1.5, 'FontWeight', 'Bold', 'FontSize', 12)

%% FRF
clear ax h
fig(3) = figure (3) ; clf
set(fig(3),'Color','w','Units','inches','Position',[4 1 3*3 3*2])
ax(1) = subplot(3,3,1); hold on ; ylabel('Gain (°/°)') ; title('Body')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,bI),...
              GRAND.fly_stats.std.IOGain.mean(:,bI), IOFreq, ...
              1, 1, bcolor, 'b', 0.2, 2);
          
ax(2) = subplot(3,3,2); hold on ; title('Head')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,hI),...
              GRAND.fly_stats.std.IOGain.mean(:,hI), IOFreq, ...
              1, 1, hcolor, 'b', 0.2, 2);
          
ax(3) = subplot(3,3,3); hold on ; title('Gaze: Body + Head')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), squeeze(GRAND.all.IOGain(:,gI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(3),h.line(3)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,gI),...
              GRAND.fly_stats.std.IOGain.mean(:,gI), IOFreq, ...
              1, 1, gcolor, 'b', 0.2, 2);       
                 
ax(4) = subplot(3,3,4); hold on ; ylabel('Phase Difference (°)')
    plot([0 10], [0 0], '--k')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,bI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(4),h.line(4)] = PlotPatch(rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,bI)),...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,bI)), IOFreq, ...
              1, 1, bcolor, 'b', 0.2, 2);
          
ax(5) = subplot(3,3,5); hold on
    plot([0 10], [0 0], '--k')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,hI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    
	phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,hI));
    %(phase>0) = phase(phase>0) - 360;
    [h.patch(5),h.line(5)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,hI)), IOFreq, ...
              1, 1, hcolor, 'b', 0.2, 2);
          
ax(6) = subplot(3,3,6); hold on
    plot([0 10], [0 0], '--k')
    plot(squeeze(GRAND.all.IOFv(:,1,:)), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,gI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    
	phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,gI));
    %(phase>0) = phase(phase>0) - 360;
    [h.patch(6),h.line(6)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,gI)), IOFreq, ...
              1, 1, gcolor, 'b', 0.2, 2);
               
ax(7) = subplot(3,3,7); hold on ; ylabel('Coherence')
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,bI),...
              GRAND.fly_stats.std.Cohr.mean(:,bI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, bcolor,'b', 0.2, 1);
%     [h.patch(7),h.line(7)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,bI),...
%               GRAND.fly_stats.std.IOCohr.mean(:,bI), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);
          
ax(8) = subplot(3,3,8); hold on
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,hI),...
              GRAND.fly_stats.std.Cohr.mean(:,hI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, hcolor,'b', 0.2, 1);
%     [h.patch(8),h.line(8)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,hI),...
%               GRAND.fly_stats.std.IOCohr.mean(:,hI), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);

ax(9) = subplot(3,3,9); hold on
    plot(squeeze(GRAND.all.Fv(:,1,:)), squeeze(GRAND.all.Cohr(:,gI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,gI),...
              GRAND.fly_stats.std.Cohr.mean(:,gI), GRAND.fly_stats.mean.Fv.mean(:,1), ...
              1, 1, gcolor, 'b', 0.2, 1);
%     [h.patch(9),h.line(9)] = PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,gI),...
%               GRAND.fly_stats.std.IOCohr.mean(:,gI), IOFreq, ...
%               1, 1, ccolor, 'b', 0.2, 2);
          
set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize', 20')
set(ax, 'LineWidth', 1.5, 'FontWeight', 'bold', 'FontSize', 12, 'XLim', [0 10],...
    'XGrid','on','YGrid','on','box','on')
set(ax(1:3),'YLim',[0 1.2])
set(ax(4:6),'YLim',[-200 100])
set(ax(7:9),'YLim',[0 1])
linkaxes(ax,'x')
linkaxes(ax(1:3),'y')
linkaxes(ax(4:6),'y')
linkaxes(ax(7:9),'y')
% set(ax(1:4),'XTickLabels',[])
YLabelHC = get(ax(7:9), 'XLabel');
set([YLabelHC{:}], 'String', 'Frequency (Hz)')
% set(ax,'XScale','log')
align_Ylabels(fig(3))

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