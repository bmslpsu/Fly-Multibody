function RMSE=CalculateRMSE_FIT(mag_out,phase_out,mag_data,phase_data,weights)
%function input:
%mag_out phase_out : simulated magnitude and phase (model)
%mage_data phase_data : experimental data
%weights : fitting weights used to estimate transfer function

%finds an error metric based on RMSE and weights of fit paramters
for i=1:length(mag_out)
    mag_outV(i)=mag_out(:,:,i);
    phase_outV(i)=phase_out(:,:,i);
end
phase_outV=phase_outV;
phase_outV=phase_outV*pi/180;

phase_data=phase_data*pi/180;
Diff_mag=sum(sqrt((mag_outV'-mag_data).^2))/length(mag_data);
Diff_phase=sum(sqrt((phase_outV'-phase_data).^2))/length(phase_data);
RMSE=Diff_mag*weights(2)+Diff_phase*(1-weights(2));