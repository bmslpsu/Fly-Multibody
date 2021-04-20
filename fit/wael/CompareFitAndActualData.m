function []  = CompareFitAndActualData(TF_intact, TF_damaged,T_intact_averaged,T_damaged_averaged,  freq)
%% load the averaged exprimental data
load('ErrorData.mat')
T_damaged_averaged=table2array(T_damaged_averaged);
T_intact_averaged=table2array(T_intact_averaged);
for i=1:length(freq)
    index_intact=find(T_intact_averaged(:,2)==freq(i));
    index_damaged=find(T_damaged_averaged(:,2)==freq(i));
    
    gain_mean_I(i)=mean(T_intact_averaged(index_intact,3));
    gain_std_I(i)=std(T_intact_averaged(index_intact,3));
    
    gain_mean_D(i)=mean(T_damaged_averaged(index_damaged,3));
    gain_std_D(i)=std(T_damaged_averaged(index_damaged,3));
    
    phase_mean_I(i)=circ_mean(T_intact_averaged(index_intact,4)*pi/180)*180/pi;
    phase_std_I(i)=circ_std(T_intact_averaged(index_intact,4)*pi/180)*180/pi;
    
    phase_mean_D(i)=circ_mean(T_damaged_averaged(index_damaged,4)*pi/180)*180/pi;
    phase_std_D(i)=circ_std(T_damaged_averaged(index_damaged,4)*pi/180)*180/pi;
    
end
%%
[gain_i, phase_i,w_out]=bode(TF_intact);
[gain_d, phase_d,w_out_damaged]=bode(TF_damaged);
for i=1:length(gain_i)
    gain_i1(i)= gain_i(1,1,i);
    phase_i1(i)= phase_i(1,1,i);
    
end

for i=1:length(gain_d)
    gain_d1(i)= gain_d(1,1,i);
    phase_d1(i)= phase_d(1,1,i);
end

subplot(2,1,1)
hold on
PlotInterval_boxplot(freq,gain_mean_I, gain_std_I, 1,'b') % plots actual data
plot(w_out/(2*pi),gain_i1,'b--') %plots simulated data from fit TF
PlotInterval_boxplot(freq,gain_mean_D, gain_std_D, 1,'r')
plot(w_out_damaged/(2*pi),gain_d1,'r--')
xlim([0.7 2.4])
legend('Intact Actual','Intact simulated','Damaged Actual','Damaged simulated')
ylabel('Gain')
subplot(2,1,2)
hold on
PlotInterval_boxplot(freq,phase_mean_I, phase_std_I, 1,'b')
plot(w_out/(2*pi),phase_i1,'b--')
PlotInterval_boxplot(freq,phase_mean_D, phase_std_D, 1,'r')
plot(w_out_damaged/(2*pi),phase_d1,'r--')
xlim([0.7 2.4])
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')