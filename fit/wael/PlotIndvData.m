function []=PlotIndvData(FRF_data,Indv_Plots)
%% OUTPUT
% plots of the raw data for each fly alone
%% INPUT
% FRF_data: frequency domain data formatted according to Ben
fly_num=size(FRF_data.ref2body.fly.gain,2); %gets the size of the
figure
subplot(2,1,1)
hold on
plot(FRF_data.IOFv{1},FRF_data.ref2body.fly.gain)
subplot(2,1,2)
hold on
plot(FRF_data.IOFv{1},FRF_data.ref2body.fly.phase)
%% individual plots
if Indv_Plots==true
    f_ind = figure;
    ax = axes(f_ind);
    ax.Units = 'pixels';
    ax.Position = [95 95 325 325];
    c = uicontrol;
    c.String = 'Plot Data';
    c.Value=false;
    for i=1:fly_num
        waitforbuttonpress;
        subplot(2,1,1)
        plot(FRF_data.IOFv{1},FRF_data.ref2body.fly.gain(:,i))
        subplot(2,1,2)
        plot(FRF_data.IOFv{1},FRF_data.ref2body.fly.phase(:,i))
        disp(['Fly number: ' num2str(i)]); 
        c.Value=false;
    end
end