function [num_bestFit,den_bestFit, delay_bestFit,RMSE_min,FitPrecent,exitFlag_bestFit,Index_rejected]...
    =FitBodePlotInformationV4_ConWithVariableDelay(npoles, nzeros,weights,delay_time,freq_FitRange,FRF_data,SpecialFit)
%V4: A good working code with exception to perform fitting for two special
%functions; leaky integrator and for a TF with a pole at 0
%% FUNCTION Description
%% INPUTS
%npoles: Number of poles
%nzeros: Number of zeros
%weights: [w1 w2] Fitting weights that define how TF is fit. w1 can be any
%real number positive values penalize errors in fitting high freq
%componenets. w2: [0 1] higher values favor phase data for fit (0.5 equal
%fit).
%delay_time: vector with range of delays to test
% freq_FitRange: frequency range of the data [freq_min freq_max]
% FRF_data: freq response data which contians IO freq, mag, and phase diff
%% SpecialFit: special fit options for certain transfer function forms
    %0 not change (standard transfer function with real number coeffs)
    %1 eliminate smallest den coeff (force one pole to be zero)
    %2 eliminate smallest two orders of den
    %3 leaky integrator

%% OUTPUS------------------------------------------------------------------
%num_bestFit: Matrix that has the numerator of the best fit TF of each fly
%den_bestFit: Matrix that has the denomenator of the best fit TF of each fly
%delay_bestFit: Vector containing the delay estimaed for each fly
%RMSE_min: Vector containing the RMSE of the each individual best fit
%FitPrecent: How much the actual data matches the fit (based on MATLABS
%tfest's definition of goodness of fit)
%exitFlag: why fmincon stopped the optimization process
%Index_rejected: Fits that were rejected due to bad 
if isempty(SpecialFit)
    SpecialFit=0;
end
%% set up the frequency domain estimate of the transfer function
syms omega
clear zeros
global a_vec b_vec %b_vec: den coeff a_vec:num coeff
a_vec = sym('a',[1 npoles+1]);
b_vec = sym('b',[1 nzeros+1]);
a_vec=fliplr(a_vec);
b_vec=fliplr(b_vec);
assume(a_vec,'real')
assume(b_vec,'real')
assume(omega,'real') %to avoid getting complex numbers when solving
eps=10^-1; % small number used here and there for initializing data
%% optimization option for Fminsearc
option=optimset('TolFun',10^-5);
option.Display='none'; %no output text to keep workspace from getting cluttered
option.MaxIter=2000;
option.MaxFunEvals=4500;
option.UseParallel=true; %parallel toolbox
%% fit the transfer function to the data and import data
allData=figure; %figures to show actual vs fit data
f_best=figure;
for j=1:size(FRF_data.err2body.fly.gain,2) %loops through each fly
    disp(['Fly Number:' num2str(j)]);
    V_er=[];
    V_array=[];
    den_values=[];
    num_values=[];
    V_val=[];
    exitFlag=[];
    RMSE=[];
    freq=[0.35 0.55 0.9 1.45 2.25 3.45 5.45 8.55 13.7]; %consider removing later on
    
    mag_ys=FRF_data.err2body.fly.gain(:,j); %extract mag and phase from exp data
    phs_ys=FRF_data.err2body.fly.phase(:,j);
    %% remove the nan values (not used anymore but a good safety to have)
    nan_index=isnan(mag_ys);
    mag_ys(nan_index)=[];
    phs_ys(nan_index)=[];
    freq(nan_index)=[];
    %% data preprocessing
    mag_ys_normalized=mag_ys; % no need to normalize since this is gain
    data = frd(mag_ys.*exp(1j*phs_ys*pi/180),freq*2*pi,1/100); %put mag and phase into data structure
    delay_est=delayest(data)/100; %built in matlab function that estimates delay
    if delay_est>0.1 %check in case  delay is too large
        disp('estimated delay too large. Delay range reset to 20-100 ms')
        delay_est=0.1;
    end
    delay_time=delay_est-0.01:0.001:delay_est+0.01;
    n_1=length(delay_time);
    i=1;
    stopFlag=[];
    while i<=n_1
        disp(['Delay:' num2str(delay_time(i))])
        %% cost function and transfer func
        [Fs, n_den, n_num]=GenerateTransferFunctionV3(npoles,nzeros,delay_time(i),SpecialFit);%transfer function is written in magnitude of it's complex component
        V_er=GenerateV(freq*2*pi, mag_ys_normalized, phs_ys*pi/180, weights,Fs);
        V_er=vpa(V_er);
        ft_v=matlabFunction(V_er);
        theta_init=eps*ones(1,n_den+n_num); %initial optiization conditions Changed with transfer function order
        %% constraints
        lower_bound=eps*ones(1,n_den+n_num)';
        upper_bound=Inf*ones(1,n_den+n_num)';
        %% begin the fit
        % %             ga(@(prms) ft_v(prms(1),prms(2),prms(3),prms(4)),4
        [theta_est_freq, V_eval,exitFlag(i), outputIteration] = ...
            fmincon(@(prms) ft_v(prms(1),prms(2),prms(3)),theta_init,...
            [],[],[],[],lower_bound,upper_bound,[],option);
        %% extract coeff
        theta_est_freq_normalized=theta_est_freq/theta_est_freq(1);
        %% TF and RMSE
        if SpecialFit==0
            den_values(i,:)=theta_est_freq_normalized(1:n_den);
            num_values(i,:)=theta_est_freq_normalized(1+n_den:end);
            tf_t2=tf(num_values(i,:),den_values(i,:),'InputDelay',delay_time(i));
        elseif SpecialFit==1
            den_values(i,:)=theta_est_freq_normalized(1:n_den);
            num_values(i,:)=theta_est_freq_normalized(1+n_den:end);
            tf_t2=tf(num_values(i,:),[den_values(i,:) 0],'InputDelay',delay_time(i));
        elseif SpecialFit==3
            tf_t2=tf([theta_est_freq(3) theta_est_freq(3)*theta_est_freq(5)+theta_est_freq(4)]...
                ,[theta_est_freq(1) theta_est_freq(1)*theta_est_freq(5)+theta_est_freq(2) ...
                theta_est_freq(2)*theta_est_freq(5)],'InputDelay',delay_time(i));
            den_values(i,:)=[theta_est_freq(1) theta_est_freq(1)*theta_est_freq(5)+theta_est_freq(2) ...
                theta_est_freq(2)*theta_est_freq(5)];
            num_values(i,:)=[theta_est_freq(3) theta_est_freq(3)*theta_est_freq(5)+theta_est_freq(4)];
        end
        [mag_out, phase_out, ~]=bode(tf_t2,freq*2*pi); %convert freq to rad
        RMSE(i)=CalculateRMSE_FIT(mag_out,phase_out,mag_ys_normalized,phs_ys,weights);
        
        i=i+1;
    end
    [RMSE_min(j), index_bestFit]=min(RMSE);
    num_bestFit(j,:)=num_values(index_bestFit,:);
    den_bestFit(j,:)=den_values(index_bestFit,:);
    delay_bestFit(j)=delay_time(index_bestFit);
    exitFlag_bestFit(j)=exitFlag(index_bestFit);
    TF_bestFit=tf(num_bestFit(j,:),den_bestFit(j,:),'InputDelay',delay_bestFit(j));
    if SpecialFit==1
        TF_bestFit=tf(num_bestFit(j,:),[den_bestFit(j,:) 0],'InputDelay',delay_bestFit(j));
    end
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
    plot(freq,mag_ys_normalized,'r')
    hold on
    plot(freq_bestFit/2/pi,mag_best_vec,'k')
    ylabel('Gain')
    subplot(2,1,2)
    plot(freq,phs_ys,'r')
    hold on
    plot(freq_bestFit/2/pi,phase_best_vec,'k')
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
index_remove=find(FitPrecent<0.4); % choice of cutoff is arbitrary
RMSE_min(index_remove)=[];
num_bestFit(index_remove,:)=[];
den_bestFit(index_remove,:)=[];
delay_bestFit(index_remove)=[];
FitPrecent(index_remove)=[];
exitFlag_bestFit(index_remove)=[];
Index_rejected=index_remove; %determine the index of the flies that were removed