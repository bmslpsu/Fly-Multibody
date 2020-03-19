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
root.benifly = fullfile(root.reg ,'tracked_headwing');
root.head = fullfile(root.reg ,'tracked_head');

% Select files
[FILES, ~] = uigetfile({'*.mat'}, 'Select trials', root.body, 'MultiSelect','on');
FILES = cellstr(FILES)';

[D,I,N,U,T,~,~,basename] = GetFileData(FILES,'',false);

%% Get Data %%
IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
Fs = 100;
Fc = 15;
tintrp = (0:(1/Fs):20)';
function_length = 20;
[b,a] = butter(2, Fc/(Fs/2),'low');
for kk = 1:N{1,end}
    disp(kk)
    % Load DAQ, body, head, & wing data
	data.daq     = load(fullfile(root.daq,  FILES{kk}),'data','t_p'); % load camera trigger & pattern x-position
    data.body    = load(fullfile(root.body, FILES{kk}),'bAngles'); % load body angles
	data.head    = load(fullfile(root.head, FILES{kk}),'hAngles'); % load head angles
 	data.benifly = ImportBenifly(fullfile(root.benifly, ...
                            [basename{kk} '.csv'])); % load head & wing angles from Benifly
    
 	% Filter wing angles
    lwing = hampel(data.benifly.Time, data.benifly.LWing);
    rwing = hampel(data.benifly.Time, data.benifly.RWing);
	lwing = rad2deg(filtfilt(b,a,lwing));
    rwing = rad2deg(filtfilt(b,a,rwing));
    
    % Sync all the signals
    daq_time = data.daq.t_p;
    daq_pattern = data.daq.data(:,2);
    trigger = data.daq.data(:,1);
    
    [TRIG,PAT] = sync_pattern_trigger(daq_time, daq_pattern, function_length, trigger, tintrp, false);
    trig_time = TRIG.time_sync;
    
    % Interpolate to so all signals have the same times
    pat = 3.75*PAT.pos;
    pat = pat - mean(pat);
    body = -data.body.bAngles;
    body = body - mean(body);
	head = -data.head.hAngles;
    head = head - mean(head);
    
    Pattern = interp1(PAT.time_sync, pat, tintrp, 'nearest');
    Body    = interp1(trig_time, body,  tintrp, 'pchip');
    Head    = interp1(trig_time, head,  tintrp, 'pchip');
    LWing   = interp1(trig_time, lwing, tintrp, 'pchip');
    RWing   = interp1(trig_time, rwing, tintrp, 'pchip');
    dWBA    = interp1(trig_time, lwing-rwing, tintrp, 'pchip');
   
    SYS = bode_sysID(tintrp, Pattern, IOFreq, true, Body, Head);
    
    pause
    close all
    
%     [Fv, Mag1 , Phs1 , FREQ1] = FFT(tintrp, Body);
%     [Fv, Mag2 , Phs2 , FREQ2] = FFT(tintrp, RWing);
%     [Fv, Mag3 , Phs3 , FREQ3] = FFT(tintrp, dWBA);
    
%     figure (1) ; clf ; hold on
%     plot(Fv, Mag3,'k','LineWidth',1)
%     plot(Fv, Mag2,'r','LineWidth',1)
%     plot(Fv, Mag1,'b','LineWidth',1)
%     xlim([0.5 10])
%     ylim([0 20])
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Magno_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end