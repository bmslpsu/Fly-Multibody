function [num_bestFit,den_bestFit, delay_bestFit,RMSE_min,FitPrecent]...
    =FitBodePlotInformationV2_BensData(npoles, nzeros,weights,delay_time,freq_FitRange,FRF_data,U)
%V2 is used to fit the error of the open loop function
%% set up the frequency domain estimate of the transfer function
syms omega
clear zeros
global a_vec b_vec
a_vec = sym('a',[1 npoles+1]);
b_vec = sym('b',[1 nzeros+1]);
assume(a_vec,'real')
assume(b_vec,'real')
assume(omega,'real')
eps=10^-5; % small number used here and there
%% optimization option for Fminsearc
option=optimset('TolFun',10^-5);
SpecialFit=0;
%% fit the transfer function to the data and import data
allData=figure;
f_best=figure;
for j=1:size(FRF_data.err2body.fly.gain,2)
    V_er=[];
    V_array=[];
    den_values=[];
    num_values=[];
    V_val=[];
    exitFlag=[];
    RMSE=[];
    freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7];

    mag_ys=FRF_data.err2body.fly.gain(:,j);
    phs_ys=FRF_data.err2body.fly.phase(:,j);
    %% remove the nan values
    nan_index=isnan(mag_ys);
    mag_ys(nan_index)=[];
    phs_ys(nan_index)=[];
    freq(nan_index)=[];
    mag_ys_normalized=mag_ys; % no need to normalize since this is gain
    data = frd(mag_ys.*exp(1j*phs_ys*pi/180),freq*2*pi,1/100);
    delay_est=delayest(data)/100;
    if delay_est>0.1 %check in case  delay is too large
        disp('estimated delay too large. Delay range reset to 20-100 ms')
        delay_est=0.1;
    end
    delay_time=delay_time(1):0.005:delay_est+0.01;
    n_1=length(delay_time);
    i=1;
    stopFlag=0;
    while i<=n_1
        disp(['Delay:' num2str(delay_time(i))])
        %%
        [Fs, n_den, n_num]=GenerateTransferFunctionV3(npoles,nzeros,delay_time(i),SpecialFit);%transfer function is written in magnitude of it's complex component
        V_er=GenerateV(freq*2*pi, mag_ys_normalized, phs_ys*pi/180, weights,Fs);
        V_er=vpa(V_er);
        ft_v=matlabFunction(V_er);
        theta_init=eps*ones(1,n_den+n_num); %initial optiization conditions Changed with transfer function order
        [theta_est_freq, V_eval,exitFlag, outputIteration] = ...
            fminsearch(@(prms) ft_v(prms(1),prms(2),prms(3),prms(4),prms(5)),theta_init,option);
        theta_est_freq_normalized=theta_est_freq/theta_est_freq(1);
        den_values(i,:)=theta_est_freq_normalized(1:n_den);
        num_values(i,:)=theta_est_freq_normalized(1+n_den:end);
%         den_values(i,:)=(theta_est_freq_normalized(1:npoles)); %for
%         special cases
%         num_values(i,:)=(theta_est_freq_normalized(npoles+1:end));
        tf_t2=tf(num_values(i,:),den_values(i,:),'InputDelay',delay_time(i));
        [mag_out, phase_out, ~]=bode(tf_t2,freq*2*pi); %convert freq to rad
        RMSE(i)=CalculateRMSE_FIT(mag_out,phase_out,mag_ys_normalized,phs_ys,weights);
        
        %                 figure(allData)
        %         bode(tf_t2)
        %         hold on
        %         xlim([0.7 2.4]*2*pi)
        if RMSE(i)>0.8 && stopFlag<=4
            i=i;
            stopFlag=stopFlag+1;
        elseif stopFlag>4
            i=i+1;
            stopFlag=0;
        else
            i=i+1;
            stopFlag=0;
        end
    end
    [RMSE_min(j), index_bestFit]=min(RMSE);
    num_bestFit(j,:)=num_values(index_bestFit,:);
    den_bestFit(j,:)=den_values(index_bestFit,:);
    delay_bestFit(j)=delay_time(index_bestFit);
    TF_bestFit=tf(num_bestFit(j,:),den_bestFit(j,:),'InputDelay',delay_bestFit(j));
    figure(f_best)
    hold on
    bode(TF_bestFit)
    [mag_best, phase_best,freq_bestFit]=bode(TF_bestFit,[freq_FitRange(1):0.1:freq_FitRange(2)]*2*pi);
    for k=1:length(mag_best)
        mag_best_vec(k)=mag_best(:,:,k);
        phase_best_vec(k)=phase_best(:,:,k);
    end
    xlim(freq_FitRange*2*pi)
    figure(allData)
    subplot(2,1,1)
    plot(freq*2*pi,mag_ys_normalized,'r')
    hold on
    plot(freq_bestFit,mag_best_vec,'k')
    ylabel('Gain')
    subplot(2,1,2)
    plot(freq*2*pi,phs_ys,'r')
    hold on
    plot(freq_bestFit,phase_best_vec,'k')
    ylabel('Phase difference')
    drawnow
    %% goodness of fit estimate
    [mag_best, phase_best,~]=bode(TF_bestFit,freq*2*pi);
    mag_best_vec_Test=[];
    phase_best_vec_Test=[]; %clear those for future fits
    for k=1:length(mag_best)
        mag_best_vec_Test(k)=mag_best(:,:,k);
        phase_best_vec_Test(k)=phase_best(:,:,k);
    end
    NMSE_mag(j)=goodnessOfFit(mag_best_vec_Test',mag_ys_normalized,'NMSE');
    NMSE_phase(j)=goodnessOfFit(phase_best_vec_Test',phs_ys,'NMSE');
    FitPrecent(j)=weights(2)*NMSE_mag(j)+(1-weights(2))*NMSE_phase(j);
    
    
end

% filter out bad fits
index_remove=find(RMSE_min>8); % choice of cutoff is arbitrary
RMSE_min(index_remove)=[];
num_bestFit(index_remove,:)=[];
den_bestFit(index_remove,:)=[];
delay_bestFit(index_remove)=[];
FitPrecent(index_remove)=[];