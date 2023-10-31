clear all
close all
clc


    %defining signal variables

T = 3;
A = 1;

T_max = 10*T;

f_s = 1000/T;
T_s = 1/f_s;

t = 0:T_s:T_max-T_s;            %time axis
N_total = length(t);            %number of samples 

N_pulse = T/T_s;                %number of samples of the pulse

s = zeros(1,N_total);           %signal 1: 1 row, N_total columns of zeroes
s(1,1:N_pulse) = A;             %signal 2: rectangular signal with amplitude A from column 1 to N_pulse


    %energy

e_analitical = A^2*T;           %analytical formula
e_time = sum(s.^2)*T_s;         %using the time axis samples


    %plotting the rectangular pulse

figure(1);

plot(t,s,'linewidth',1);
grid on;
xlabel('time');
ylabel('s(t)');
title(sprintf('Rectangular pulse \nA = %d   T = %d \nPulse energy = %.3f   \nEnergy on time axis = %.3f\n', A, T, e_analitical, e_time));


    %fft and spectrum

f_res = f_s/N_total;            %frequency resolution, for no interference
f = 0:f_res:f_s-f_res;          %frequency axis
f_symm = f-f_s/2;               %make the frequency axis symmetric

S = fft(s)*T_s;                 %fft of s scaled by T_s
M = abs(S);                     %magnitude of S
M_shifted = fftshift(M);        %centered wrt 0 -> symmetric

e_frequency = sum(M.^2)*f_res;


    %plotting the rectangular pulse frequency spectrum

figure(2);

plot(f_symm,M_shifted,'linewidth',1);
grid on;
xlabel('frequency');
ylabel('|S(f)|^2');
xlim([-15 15]);
title(sprintf('Rectangular pulse \nA = %d   T = %d \nPulse energy = %.3f   \nEnergy on frequency axis = %.3f\n', A, T, e_analitical, e_frequency));


    %energy percentage in first 10 lobes

e_10 = zeros(1,10);
e_10(1) = (M(1).^2+2*sum(M(2:(1/T)/f_res).^2))*f_res;

for x = 2:10
    
    e_10(x) = e_10(x-1) + 2*sum((M(((x-1)/T/f_res)+1:x/T/f_res)).^2)*f_res;
    
end

e_10_percentage = e_10./e_frequency.*100;


    %energy percentage in first 100 lobes
    
e_100 = zeros(1,100);
e_100(1:10) = e_10;

for z = 11:100
    
    e_100(z) = e_100(z-1) + 2*sum((M(((z-1)/T/f_res)+1:z/T/f_res)).^2)*f_res;
    
end

e_100_percentage = e_100./e_frequency.*100;


    %plotting energy in first 10 and 100 lobes

figure(3);

xaxis = 1:1:100;
subplot(2,1,1);
plot(xaxis(1:10),e_10_percentage,'-o');
grid on;
xlabel('number of lobes');
ylabel('energy percentage');
xlim([0 10]);
title(sprintf('A = %d \nT = %d \n\nFirst 10 lobes:', A, T));

subplot(2,1,2);
plot(xaxis,e_100_percentage,'linewidth',1);
grid on;
xlabel('number of lobes');
ylabel('energy percentage');
xlim([0 100]);
title(sprintf('First 100 lobes:'));
