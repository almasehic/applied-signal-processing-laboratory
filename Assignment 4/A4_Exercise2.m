clear all 
close all
clc


    %defining variables
    
A1 = 1;
f1 = 1;
A2 = 2;
f2 = 4;

f_s = 100;
T_s = 1/f_s;
T = 100;

t = 0:T_s:T-T_s;
N = length(t);
f_res = f_s/N;
f = 0:f_res:f_s-f_res;
f_symm = f-f_s/2;

phi1 = rand()*2*pi;
s1 = A1*cos(2*pi*f1*t+phi1);

phi2 = rand()*2*pi;
s2 = A2*cos(2*pi*f2*t+phi2);

s = s1+s2;

S = fft(s)*T_s;             %fourier transform of s
M = abs(S);                 %magnitude of S
M_shift = fftshift(M);


    %plotting
    
figure
subplot(2,1,1);
plot(t, s);
title('s(t)');
xlabel('time t');
xlim([0 100]);
ylim([-3 3]);
subplot(2,1,2);
plot(f_symm, M_shift);
title('abs(S(f))');
xlabel('frequency f');
xlim([-50 50]);


    %new signal s'
    
s_prime = rectangularPulse(0,T/2,t).*s1 + rectangularPulse(T/2,T,t).*s2;
S_Prime = fft(s_prime)*T_s;
MP = abs(S_Prime);
MP_shift = fftshift(MP);

figure
subplot(2,1,1);
plot(t, s_prime);
title('s''(t)');
xlabel('time t');
xlim([0 100]);
ylim([-2 2]);
subplot(2,1,2);
plot(f_symm, MP_shift);
title('abs(S''(f))');
xlabel('frequency f');
xlim([-50 50]);
ylim([0 50]);


    %plotting the spectrograms
    
N_fft = N/2;
N_overlap = 0;
window = ones(1,N_fft);
[S2,f2,t2] = spectrogram(s,window,N_overlap,N_fft,f_s,'centered','yaxis');

figure
imagesc(t2,f2,abs(S2)*T_s);
axis xy;
ylim([0 10]);
h = colorbar;
ylabel(h , 'abs(S)');
xlabel('time [s]');
ylabel('frequency [Hz]');

[Sp,fp,tp] = spectrogram(s_prime,window,N_overlap,N_fft,f_s,'centered','yaxis');

figure
imagesc(tp,fp,abs(Sp)*T_s);
axis xy;
ylim([0 10]);
h = colorbar;
ylabel(h , 'abs(S'')');
xlabel('time [s]');
ylabel('frequency [Hz]');


%work done by: Alma Sehic
%              s274208

