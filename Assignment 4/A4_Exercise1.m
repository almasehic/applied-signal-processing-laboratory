clear all
close all
clc


    %recording audio and defining variables
    
f_s = 44100;
T_s = 1/f_s;
T = 3;
myVoice = audiorecorder(f_s,24,1);
disp('Start recording.....');
recordblocking(myVoice, T);
disp('End of recording. Playing back...');
play(myVoice);
y = getaudiodata(myVoice);

t = 0:T_s:T-T_s;
N = length(t);
f_res = f_s/N;
f = 0:f_res:f_s-f_res;
f_symm = f-f_s/2;

Y = fft(y)*T_s;
M = abs(Y);
M_shift = fftshift(M);
max_Y = max(M_shift);
n_Y_shift = M_shift/max_Y;


    %plotting
    
figure
subplot(1,2,1);
plot(t, y);
grid on
ylabel('s(t)');
xlabel('time [s]');

subplot(1,2,2);
plot(f_symm, n_Y_shift);
grid on
ylabel('|S(f)|');
xlabel('frequency [Hz]');


%work done by: Alma Sehic
%              s274208

