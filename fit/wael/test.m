%Code now works. Modified for Ben's data which has a different formate
%compared to mine. Has fmincon now
%% this is a working code that can fit transfer function to frequency domain data
clc
clear all
close all
%% Initialize the variables for data analysis
% output data format bn(s^n-1)+b(n-1)s^n-2.../an*s^n+a(n-1)*s^(n-1)...
% the first n of theta are the the den values which start from the last one
% which is the highest order value
npoles=1;
nzeros=0;
weights=[-0.1 0.5]; % fitting weights [w1 w2]. w1: -inf to inf. w2: [0 1]
delay_time=[0.018:0.001:0.08]; % a vector containing all the possible delays
freq_range=[0.35 14]; % range of plotting and fitting
%% load the data used in the analysis. It is in the form of a table( needs gain, phase, freq, and fly number)
%removed some flies that had a strange response or bad videos
load('E:\DATA\Magno_Data\Multibody\processed\wael\FRF_SOS_HeadFree_vel_52_position.mat')
Frf_D=FRF_data;
U_D=U;
Frf_D.ref2body.fly.gain(:,[10 33])=[];
Frf_D.ref2body.fly.phase(:,[10 33])=[];
Frf_D.err2body.fly.gain(:,[10 33])=[];
Frf_D.err2body.fly.phase(:,[10 33])=[];
fly_numD=U_D.fly{1};
fly_numD([10 33])=[];
U_D.fly=mat2cell(fly_numD,length(fly_numD));
load('E:\DATA\Magno_Data\Multibody\processed\wael\FRF_SOS_HeadFree_vel_52_wing_damage_position.mat')
Frf_I=FRF_data;
U_I=U;
%% seperate the bode plots by area
root_area='S:\Public\Wael\Chapter4\DamagedWingExp\SOS_new\Wing_Pictures\AnalyzedImages\';
p_values=GenerateBodePlotByArea(root_area,3,U_D,Frf_D);
%% determine which fits to reject
PlotIndvData(Frf_I,false)
PlotIndvData(Frf_D,false)
%% fit type 
%0 not change
%1 eliminate smallest den coeff
%2 eliminate smallest two orders of den (not done yet)
%3 leaky integrator
%% the main fitting function (old fitting code but useful as a reference due to good fitting flexibility)
%usually used for testing fits before moving on, though lately it doesn't
%seem to have much use
[num_bestFitI,den_bestFitI, delay_bestFitI,RMSE_minI,FitPrecentI]...
    =FitBodePlotInformationV2_BensData(npoles, nzeros,weights,delay_time,freq_range,Frf_I,U_I);

[num_bestFitD,den_bestFitD, delay_bestFitD,RMSE_minD,FitPrecentD]...
    =FitBodePlotInformationV2_BensData(npoles, nzeros,weights,delay_time,freq_range,Frf_D,U_D);
%% try fmincon fitting code mostly used for higher order fits
%latest fitting function
SpecialFit=0; %modify to use special functions
[num_bestFitI4,den_bestFitI4, delay_bestFitI4,RMSE_minI4,FitPrecentI4,exitFlag_bestFitI4,index_removeI4]=...
    FitBodePlotInformationV4_ConWithVariableDelay(npoles, nzeros,weights,delay_time,freq_range,Frf_I,SpecialFit);
[num_bestFitD4,den_bestFitD4, delay_bestFitD4,RMSE_minD4,FitPrecentD4,exitFlag_bestFitD4,index_removeD4]...
    =FitBodePlotInformationV4_ConWithVariableDelay(npoles, nzeros,weights,delay_time,freq_range,Frf_D,SpecialFit);
%% Compare fits to actual data
% add a 0 to the den vector when using special fit 1
TF_intact=tf(mean(num_bestFitI4),[mean(den_bestFitI4)],'InputDelay',mean(delay_bestFitI4));
TF_damaged=tf(mean(num_bestFitD4),[mean(den_bestFitD4)],'InputDelay',mean(delay_bestFitD4));
%% plot actual averaged vs fit
figure
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
CompareFitAndActualDataV2( TF_intact,freq,Frf_I,'b','--b')
CompareFitAndActualDataV2( TF_damaged,freq,Frf_D,'r','--r')
%% compare the gain and phase
figure
subplot(2,1,1)
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
gain_I=Frf_I.ref2body.fly.gain;
Phase_I=Frf_I.ref2body.fly.phase;
PlotInterval(freq, mean(gain_I,2)', Frf_I.ref2body.grand_std.gain', 1,'b')
subplot(2,1,2)
PlotInterval(freq,mean(Phase_I,2)', Frf_I.ref2body.grand_std.phase', 1,'b')

subplot(2,1,1)
gain_D=Frf_D.ref2body.fly.gain;
PlotInterval(freq,mean(gain_D,2)', Frf_D.ref2body.grand_std.gain', 1,'r')
ylabel('Gain')
ylim([-0 1])
subplot(2,1,2)
phase_D=Frf_D.ref2body.fly.phase;
PlotInterval(freq,mean(phase_D,2)', Frf_D.ref2body.grand_std.phase', 1,'r')
xlabel('Frequency (Hz)')
ylabel('Phase difference')
ylim([-250 25])
%% error compensation plot
figure
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
PlotInterval_boxplot(freq, Frf_I.ref2body.grand_mean.error, Frf_I.ref2body.grand_std.error, 1,'b')
PlotInterval_boxplot(freq,Frf_D.ref2body.grand_mean.error, Frf_D.ref2body.grand_std.error, 1,'r')
ylabel('Gain')
xlabel('Frequency (Hz)')
title('Compensation error')
%% coherence plots
figure
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
gain_I=Frf_I.ref2body.fly.IO_coherence;
PlotInterval(freq, mean(gain_I,2)', Frf_I.ref2body.grand_std.IO_coherence', 1,'b')
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
gain_D=Frf_D.ref2body.fly.IO_coherence;
PlotInterval(freq,mean(gain_D,2)', Frf_D.ref2body.grand_std.IO_coherence', 1,'r')
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])

%% compare the gains and phase of both flies
[p_mag, p_phase]=ConductTtestBen(Frf_I,Frf_D);
%% coherence t-test
[p_coh,f_coh]=ConductTtestBen_Coh(Frf_I,Frf_D);
%% compare controller values
[p_KP, p_Damping]=CreateControllerBoxplot(num_bestFitI4,den_bestFitI4,num_bestFitD4,den_bestFitD4);
%% generate individual bode plots
Fly_count_intact=GenerateIndvBodePlot(num_bestFitI,den_bestFitI,delay_bestFitI,freq);
suptitle('Intact-wing OL bode plots')
Fly_count_damaged=GenerateIndvBodePlot(num_bestFitD,den_bestFitD,delay_bestFitD,freq);
suptitle('Damaged-wing OL bode plots')
%% POLEZERO MAP
figure
pzmap(TF_intact)
hold on
pzmap(TF_damaged)
%% Head Analysis--------------------------------
%% conduct head analysis
load("E:\DATA\Magno_Data\Multibody\processed\FRF_SOS_HeadFree_vel_52_position.mat")
Frf_I=FRF_data;
U_I=U;
load("E:\DATA\Magno_Data\Multibody\processed\FRF_SOS_HeadFree_vel_52_wing_damage_position.mat")
Frf_D=FRF_data;
U_D=U;
Frf_D.ref2body.fly.gain(:,[10 33])=[];
Frf_D.ref2body.fly.phase(:,[10 33])=[];
Frf_D.err2body.fly.gain(:,[10 33])=[];
Frf_D.err2body.fly.phase(:,[10 33])=[];
fly_numD=U_D.fly{1};
fly_numD([10 33])=[];
U_D.fly=mat2cell(fly_numD,length(fly_numD));
%%
figure
subplot(2,1,1)
freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];
gain_I=Frf_I.ref2head.fly.gain;
Phase_I=Frf_I.ref2head.fly.phase;
PlotInterval(freq, mean(gain_I,2)', Frf_I.ref2head.grand_std.gain', 1,'b')
subplot(2,1,2)
PlotInterval(freq,mean(Phase_I,2)', Frf_I.ref2head.grand_std.phase', 1,'b')

subplot(2,1,1)
gain_D=Frf_D.ref2head.fly.gain;
PlotInterval(freq,mean(gain_D,2)', Frf_D.ref2head.grand_std.gain', 1,'r')
ylabel('Gain')
ylim([-0 1])
subplot(2,1,2)
phase_D=Frf_D.ref2head.fly.phase;
PlotInterval(freq,mean(phase_D,2)', Frf_D.ref2head.grand_std.phase', 1,'r')
xlabel('Frequency (Hz)')
ylabel('Phase difference')
ylim([-150 150])
%%
for i=1:9
    damaged_mag=Frf_D.ref2head.fly.gain(i,:);
    intact_mag=Frf_I.ref2head.fly.gain(i,:);
    [~,p_mag_Head(i)]=ttest2(intact_mag,damaged_mag);
    
    damaged_phase=Frf_D.ref2head.fly.phase(i,:);
    intact_phase=Frf_I.ref2head.fly.phase(i,:);
    [~,p_phase_Head(i)]=ttest2(intact_phase,damaged_phase);
    
end
%% FUNCTIONS-------------
function [p_KP, p_Damping]=CreateControllerBoxplot(num_bestFitI,den_bestFitI,num_bestFitD,den_bestFitD)
%% control coeff
figure
subplot(2,1,1)
boxplot([num_bestFitI(:,1); num_bestFitD(:,1)],[num_bestFitI(:,1)*0; 1+0*num_bestFitD(:,1)],'Labels',["Intact","Damaged"])
title('P')
subplot(2,1,2)
boxplot([num_bestFitI(:,2); num_bestFitD(:,2)],[num_bestFitI(:,2)*0; 1+0*num_bestFitD(:,2)],'Labels',["Intact","Damaged"])
title('I')
[~,p_KP]=ttest2(num_bestFitI(:,1), num_bestFitD(:,1));
[~,p_KI]=ttest2(num_bestFitI(:,2), num_bestFitD(:,2));
suptitle('Controller coefficients')
%% damping
figure
Inertia=5.2*10^-13;
boxplot([den_bestFitI(:,2)*Inertia; den_bestFitD(:,2)*Inertia],[den_bestFitI(:,2)*0; 1+0*den_bestFitD(:,2)],'Labels',["Intact","Damaged"])
[~,p_Damping]=ttest2(den_bestFitI(:,2), den_bestFitD(:,2));
title('Damping coefficient')
end

function []  = CompareFitAndActualDataV2(TF_damaged,freq,FRF_data,color,color2)
%% load the averaged exprimental data
[gain_d, phase_d,w_out_damaged]=bode(TF_damaged,freq*2*pi);
for i=1:length(gain_d)
    gain_d1(i)= gain_d(1,1,i);
    phase_d1(i)= phase_d(1,1,i);
end
gain_mean_D=FRF_data.err2body.grand_mean.gain;
phase_mean_D=FRF_data.err2body.grand_mean.phase;
gain_std_D=FRF_data.err2body.grand_std.gain;
phase_std_D=FRF_data.err2body.grand_std.phase;
subplot(2,1,1)
hold on
PlotInterval(freq,gain_mean_D', gain_std_D', 1,color)
plot(w_out_damaged/(2*pi),gain_d1,color2)
xlim([0 14])
ylabel('Gain')
subplot(2,1,2)
hold on
PlotInterval(freq,phase_mean_D', phase_std_D', 1,color)
plot(w_out_damaged/(2*pi),phase_d1,color2)
xlim([0 14])
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
end

function [p_mag, p_phase]=ConductTtestBen(FrF_I,Frf_D)

for i=1:9
    damaged_mag=Frf_D.ref2body.fly.gain(i,:);
    intact_mag=FrF_I.ref2body.fly.gain(i,:);
    [~,p_mag(i)]=ttest2(intact_mag,damaged_mag);
    
    damaged_phase=Frf_D.ref2body.fly.phase(i,:);
    intact_phase=FrF_I.ref2body.fly.phase(i,:);
    [~,p_phase(i)]=ttest2(intact_phase,damaged_phase);
    
end
end

function [p_coh, f_coh]=ConductTtestBen_Coh(FRF_I,FRF_D)

for i=1:9
    damaged_mag=FRF_D.ref2body.fly.IO_coherence(i,:);
    intact_mag=FRF_I.ref2body.fly.IO_coherence(i,:);
    [~,p_coh(i)]=ttest2(intact_mag,damaged_mag);
    [~,f_coh(i)]=vartest2(intact_mag,damaged_mag);
end
end