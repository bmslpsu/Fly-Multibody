function [] = MakeData_SOS_v1_HeadFree(rootdir)
%% MakeData_SOS_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
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
    
    [TRIG,PAT]  = sync_pattern_trigger(daq_time, daq_pattern, function_length, trigger, tintrp, false);
    trig_time   = TRIG.time_sync;
    
  	% Filter wing angles
    lwing = hampel(data.benifly.Time, data.benifly.LWing);
    rwing = hampel(data.benifly.Time, data.benifly.RWing);
	lwing = rad2deg(filtfilt(b,a,lwing));
    rwing = rad2deg(filtfilt(b,a,rwing));
    
    % Get pattern, head, & body anlges
    pat = 3.75*PAT.pos;
    pat = pat - mean(pat);
    body = data.body.bAngles;
    body = body - mean(body);
	head = data.head.hAngles;
    % head = rad2deg(-data.benifly.Head);
    head = head - mean(head);
    
    % Interpolate so all signals have the same times
    Pattern = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
    
    SYS_stim2_head_body = bode_sysID(tintrp, Pattern, IOFreq, debug, Body, Head);
    SYS_stim2_wing = bode_sysID(tintrp, Pattern, IOFreq, debug, LWing, RWing, dWBA);
    
    SYS_all = CatStructFields(2, SYS_stim2_head_body, SYS_stim2_wing);
    
    ALL{I.fly(kk)}(end+1,1) = SYS_stim2_head_body;
end

%% Group Data
clc
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

%% Coherence
bI = 1;
hI = 2;

clear ax
fig(1) = figure (1) ; clf
set(fig(1),'Color','w','Units','inches','Position',[2 2 3 5])
ax(1) = subplot(2,1,1); hold on ; ylabel('Body Coherence')
    plot(squeeze(GRAND.all.Fv), squeeze(GRAND.all.Cohr(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,bI),...
              GRAND.fly_stats.std.Cohr.mean(:,bI), GRAND.fly_stats.mean.Fv.mean, ...
              1, 1, 'k','b', 0.2, 2);
    PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,bI),...
              GRAND.fly_stats.std.IOCohr.mean(:,bI), IOFreq, ...
              1, 1, 'r','r', 0.2, 2);
          
ax(2) = subplot(2,1,2); hold on ; ylabel('Head Coherence') ; xlabel('Frequency (Hz)')
    plot(squeeze(GRAND.all.Fv), squeeze(GRAND.all.Cohr(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.Cohr.mean(:,hI),...
              GRAND.fly_stats.std.Cohr.mean(:,hI), GRAND.fly_stats.mean.Fv.mean, ...
              1, 1, 'k','b', 0.2, 2);
    PlotPatch(GRAND.fly_stats.mean.IOCohr.mean(:,hI),...
              GRAND.fly_stats.std.IOCohr.mean(:,hI), IOFreq, ...
              1, 1, 'r','r', 0.2, 2);

set(ax,'LineWidth',1.5,'FontWeight','Bold','FontSize',10,'XLim',[0 10],'YLim',[0 1])
linkaxes(ax)
% YLabelHC = get(ax, 'YLabel');
% set([YLabelHC{:}], 'String', 'Head Position (°)')

%% Complex Gain
clear ax
cc = hsv(length(IOFreq));
fig(2) = figure (2) ; clf
set(fig(2),'Color','w','Units','inches','Position',[2 2 8 5])
ax(1) = subplot(1,2,1); hold on ; title('Body Complex Gain')
    Real = real(squeeze(GRAND.all.IOFRF(:,bI,:)))';
    Img  = imag(squeeze(GRAND.all.IOFRF(:,bI,:)))';
    [~,~] = ComplexAxes(0:0.2:1.2);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
ax(2) = subplot(1,2,2); hold on ; title('Head Complex Gain')
    Real = real(squeeze(GRAND.all.IOFRF(:,hI,:))');
    Img  = -imag(squeeze(GRAND.all.IOFRF(:,hI,:))');
    [~,~] = ComplexAxes(0:0.2:0.5);
    hh = plot(Real,Img,'.','MarkerSize',10);
    for c = 1:length(hh)
       hh(c).Color = cc(c,:); 
    end
    legend(hh,string(IOFreq) + " Hz",'Box','off','Location','northoutside')
    
set(ax,'LineWidth',1.5,'FontWeight','Bold','FontSize',10)

%% FRF
clear ax h
fig(3) = figure (3) ; clf
set(fig(3),'Color','w','Units','inches','Position',[2 2 6 5])
ax(1) = subplot(2,2,1); hold on ; ylabel('Gain (°/°)') ; title('Body')
    plot(squeeze(GRAND.all.IOFv), squeeze(GRAND.all.IOGain(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(1),h.line(1)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,bI),...
              GRAND.fly_stats.std.IOGain.mean(:,bI), IOFreq, ...
              1, 1, 'k','b', 0.2, 2);
          
ax(2) = subplot(2,2,2); hold on ; title('Head')
    plot(squeeze(GRAND.all.IOFv), squeeze(GRAND.all.IOGain(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(2),h.line(2)] = PlotPatch(GRAND.fly_stats.mean.IOGain.mean(:,hI),...
              GRAND.fly_stats.std.IOGain.mean(:,hI), IOFreq, ...
              1, 1, 'k','b', 0.2, 2);
          
ax(3) = subplot(2,2,3); hold on ; ylabel('Phase Difference (°)')
    plot(squeeze(GRAND.all.IOFv), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,bI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    [h.patch(3),h.line(3)] = PlotPatch(rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,bI)),...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,bI)), IOFreq, ...
              1, 1, 'k','b', 0.2, 2);
ax(4) = subplot(2,2,4); hold on
    plot(squeeze(GRAND.all.IOFv), rad2deg(squeeze(GRAND.all.IOPhaseDiff(:,hI,:))), ...
        'Color', [0.5 0.5 0.5 0.3])
    
	phase = rad2deg(GRAND.fly_stats.circ_mean.IOPhaseDiff.circ_mean(:,hI));
    %(phase>0) = phase(phase>0) - 360;
    [h.patch(4),h.line(4)] = PlotPatch(phase,...
              rad2deg(GRAND.fly_stats.circ_std.IOPhaseDiff.circ_mean(:,hI)), IOFreq, ...
              1, 1, 'k','b', 0.2, 2);

set(h.line,'Marker','.','MarkerFaceColor','none','MarkerSize',25')
set(ax,'LineWidth',1.5,'FontWeight','Bold','FontSize',10,'XLim',[0 10],...
    'XGrid','on','YGrid','on')
set(ax(1:2),'YLim',[0 1.2])
set(ax(3:4),'YLim',[-300 200])
linkaxes(ax,'x')
linkaxes(ax(1:2),'y')
linkaxes(ax(3:4),'y')
YLabelHC = get(ax(3:4), 'XLabel');
set([YLabelHC{:}], 'String', 'Frequency (Hz)')
% set(ax,'XScale','log')

%% Time
clear ax
fig(4) = figure (4) ; clf
set(fig(4),'Color','w','Units','inches','Position',[2 2 10 6])
ax(1) = subplot(2,1,1); hold on ; ylabel('Body (°)')
    plot(squeeze(GRAND.all.Time), squeeze(GRAND.all.State(:,bI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.State.mean(:,bI),...
              GRAND.fly_stats.std.State.mean(:,bI), GRAND.fly_stats.mean.Time.mean, ...
              1, 1, 'k','b', 0.2, 2);
          
ax(2) = subplot(2,1,2); hold on ; ylabel('Head (°)') ; xlabel('Time (s)')
    plot(squeeze(GRAND.all.Time), squeeze(GRAND.all.State(:,hI,:)), ...
        'Color', [0.5 0.5 0.5 0.3])
  	PlotPatch(GRAND.fly_stats.mean.State.mean(:,hI),...
              GRAND.fly_stats.std.State.mean(:,hI), GRAND.fly_stats.mean.Time.mean, ...
              1, 1, 'k','b', 0.2, 2);

set(ax,'LineWidth',1.5,'FontWeight','Bold','FontSize',10,'XLim',[0 20])
linkaxes(ax,'x')


%% SAVE %%
disp('Saving...')
save(['H:\DATA\Magno_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'], ...
    'FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end