function [Fs,n_den, n_num]=GenerateTransferFunctionV3(poles,zeros,delay_time,specialFunction)
global a_vec b_vec tau
num=0;
den=0;
syms omega


for j=0:poles
    den=den+a_vec(j+1)*(1i*omega)^(j);
end

for j=0:zeros
    num=num+b_vec(j+1)*(1i*omega)^j;
end

Fs=num/den*exp(-delay_time*omega*1i);
%% special function section
if specialFunction==0
    n_den=length(a_vec);
    n_num=length(b_vec);
elseif specialFunction==1
    Fs=subs(Fs,a_vec(1),0); % removes the lowest coeff denominator to test certain types of fitting
    n_den=length(a_vec)-1;
    n_num=length(b_vec);
elseif specialFunction==2
    Fs=subs(Fs,a_vec(1),0);
    Fs=subs(Fs,a_vec(2),0);
        n_den=length(a_vec)-2;
    n_num=length(b_vec);
elseif specialFunction==3
    b_vec = sym('b',[1 3]);
    Fs=exp(-delay_time*omega*1i)*(b_vec(1)+b_vec(2)/(omega*1i+b_vec(3)))...
        /(a_vec(3)*omega*1i+a_vec(2));
    n_den=2;
    n_num=3;
end