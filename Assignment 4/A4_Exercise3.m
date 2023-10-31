clear all
close all
clc


    %imporing audiofile and defining variables
    
filename = 'ex5_5_A_major_scale.wav';    
[s, f_s] = audioread(filename);
playObj = audioplayer(s,f_s)
playblocking(playObj)

N = length(s);
T_s = 1/f_s;
T = N/f_s;
t = 0:T_s:T-T_s;

f_res = f_s/N;
f = 0:f_res:f_s-f_res;
f_symm = f-f_s/2;

S = fft(s)*T_s;
M = abs(S);
M_s = fftshift(M);
    

    %plotting
    
figure
plot(t, s);
grid on
title('A major scale');
xlim([0 10]);
xlabel('time [s]');
ylabel('amplitude');

figure 
plot(f_symm, M_s);
title('abs(S)');
xlim([-800 800]);
xlabel('ff');


    %spectrogram
    
N_fft = 2^14;
N_overlap = 8000;
window = hamming(10000);

[S1,f1,t1] = spectrogram(s,window,N_overlap,N_fft,f_s,'centered','yaxis');

figure
imagesc(t1,f1,abs(S1)*T_s);
colormap('turbo')
axis xy;
ylim([0 800]);
h = colorbar;
ylabel(h , 'abs(S)');
xlabel('time [s]');
ylabel('frequency [Hz]');

s2 = shiftPitch(s,3);

[S2,f2,t2] = spectrogram(s2,window,N_overlap,N_fft,f_s,'centered','yaxis');

figure
imagesc(t2,f2,abs(S2)*T_s);
colormap('turbo')
axis xy;
ylim([0 800]);
h = colorbar;
ylabel(h , 'abs(S2)');
xlabel('time [s]');
ylabel('frequency [Hz]');


%work done by: Alma Sehic
%              s274208

