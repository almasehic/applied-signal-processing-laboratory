clear all
close all
clc


    %defining variables
  
T_max = 1;
f_s = 100;                  %sampling frequency
T_s = 1/f_s;                %sampling period
t = 0:T_s:T_max-T_s;        %time axis
N_total = length(t);        %number of samples
f_res = f_s/N_total;
f_0 = 1;
A = 1;

s = A*cos(2*pi*f_0*t);      %signal s(t);

m = 3;                                      %quantization bits
delta = (max(s)-min(s))/(2^m);              %average value
partition = [-((2^m)/2-1)*delta:delta:((2^m)/2-1)*delta];
codebook = [-((2^m)/2-1)*delta-delta/2:delta:((2^m)/2-1)*delta+delta/2];
[index,sq] = quantiz(s,partition,codebook);

    
    %first plot

figure(1)
plot(t,s);
hold on
plot(t,sq);
legend('s(t)','sq(t)')
grid on
title(sprintf('original signal s(t) vs. quantized signal s_q(t)\n\nm=%d    A=%d   f0=%d',m,A,f_0))

S = fft(s)*T_s;             %fourier transform of the signal
SM = abs(S);                %magnitude of S(f)
M = 20*log10(SM);           %dB
S_shift = fftshift(M);      %shifting for symmetry      

SQ = fft(sq)*T_s;           %fourier transform of the quantized signal
SMQ = abs(SQ);              %magnitude of Sq(f)
MQ = 20*log10(SMQ);         %dB
SQ_shift = fftshift(MQ);    %shifting for symmetry


E_s = max(SMQ);              %energy of the frequency component f_0
New_s = SM - SMQ;
E_n = sum(New_s.^2)*f_res;   %energy of the noise

SNR1 = 10*log10(E_s/E_n);         %calculating the SNR using the energy ratio
SNR2 = 10*log10(1.5*4^m);         %calculating the SNR using a derived formula 

f_res = f_s/N_total;            %frequency resolution
f = 0:f_res:f_s-f_res;          %frequency axis
f_symm = f-f_s/2;               %symmetric axis


    %second plot

figure(2)
subplot(2,1,1);
plot(f_symm,S_shift);
grid on
xlabel('f');
ylabel('|S(f)|^2[dB]');
xlim([-50 50])
ylim([-50 0])
grid on

subplot(2,1,2);
plot(f_symm,SQ_shift);
xlabel('f');
ylabel('|S_q(f)|^2[dB]');
xlim([-50 50])
ylim([-50 0])
grid on
title(sprintf('m = %d   SNR = %.2f [dB]', m, SNR1));



