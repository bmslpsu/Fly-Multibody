function V_er=GenerateV(freq, gain, phase, weights,Fs)
%% creates the error which will be used to find the components
V_er=0;
syms omega
digits(5)
for j=1:length(freq)
    try
        Fs_phase=angle(subs(Fs,omega,freq(j))); %phase angle of the fit transfer function
        Fs_mag=norm(subs(Fs,omega,freq(j))); %magnitude of the fit transfer function
        part1=log(Fs_mag/gain(j))^2*weights(2);  %F_fit(mag)/F_exp(mag)
        part2=(1-weights(2))*(Fs_phase-phase(j))^2; %angle(F_fit(mag)/F_exp(mag)
        V_er=V_er+(freq(j)^weights(1))*(part1+part2);
    catch % this should happen much but some trials don't have all the frequency data 
        %especially at the end
        disp('Freq data doesnt exist for this trial')
    end
end
