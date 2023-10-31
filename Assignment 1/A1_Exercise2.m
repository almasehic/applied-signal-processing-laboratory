clear all 
close all
clc


  %variable and signal definitions

T_max = 10;
f_s = 100;                      %sampling frequency
T_s = 1/f_s;                    %sampling period

fA = 1;                         %frequencies for a b c
fB = 2;
fC = 3;

phiA = rand()*2*pi;             %random phases for a b c
phiB = rand()*2*pi;
phiC = rand()*2*pi;

t = 0:T_s:T_max-T_s;            %time axis
N_total = length(t);            %total number of samples in 10 seconds

sA = sin(2*pi*fA*t + phiA);     %sine signals a b c
sB = sin(2*pi*fB*t + phiB);
sC = sin(2*pi*fC*t + phiC);

x = sA + sB + sC;                   %sum of sine signals

f_res = f_s/N_total;            %frequency resolution
f = 0:f_res:f_s-f_res;          %frequency axis
f_symm = f-f_s/2;               %symmetric frequency axis

X = fft(x)*T_s;             %fourier transform of x
M = abs(X);                 %magnitude of X
M_shift = fftshift(M);      %shifting on the f axis


  %plot

figure(1);
subplot(2,1,1);
plot(t,x,'linewidth',1);
grid on;
ylim([-3 3.5]);
xlim([0 5]);
xlabel('time t');
ylabel('x(t)');


subplot(2,1,2);
plot(f_symm,M_shift,'linewidth',1);
grid on;
xlim([-10 10]);
xlabel('frequency f');
ylabel('|X(f)|');


  %low-pass filter design

H_A = rectangularPulse(-1.4,1.4,f_symm);
Y_A = fftshift(X).*H_A;
M_A = abs(Y_A);

figure(2)
subplot(2,2,1);
plot(f_symm,H_A);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|H_A(f)|');
subplot(2,2,2);
plot(f_symm,M_A);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|Y_A(f)|');

Y_A_shift = ifftshift(Y_A);
y_a = ifft(Y_A_shift)*f_s;

subplot(2,2,3);
plot(t,y_a);
xlim([0 4]);
xlabel('time t');
ylabel('y_A(t)');
subplot(2,2,4);
plot(t,sA);
xlim([0 4]);
xlabel('time t');
ylabel('s_A(t)');


    %band-pass filter design
    
H_B = rectangularPulse(-2.4,-1.6,f_symm)+rectangularPulse(1.6,2.4,f_symm);
Y_B = fftshift(X).*H_B;
M_B = abs(Y_B);

figure(3)
subplot(2,2,1);
plot(f_symm,H_B);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|H_B(f)|');
subplot(2,2,2);
plot(f_symm,M_B);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|Y_B(f)|');

Y_B_shift = ifftshift(Y_B);
y_b = ifft(Y_B_shift)*f_s;

subplot(2,2,3);
plot(t,y_b);
xlim([0 4]);
xlabel('time t');
ylabel('y_B(t)');
subplot(2,2,4);
plot(t,sB);
xlim([0 4]);
xlabel('time t');
ylabel('s_B(t)');


    %high-pass filter design
    
H_C = 1-rectangularPulse(-2.6,2.6,f_symm);
Y_C = fftshift(X).*H_C;
M_C = abs(Y_C);

figure(4)
subplot(2,2,1);
plot(f_symm,H_C);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|H_C(f)|');
subplot(2,2,2);
plot(f_symm,M_C);
xlim([-4 4]);
xlabel('frequency f');
ylabel('|Y_C(f)|');

Y_C_shift = ifftshift(Y_C);
y_c = ifft(Y_C_shift)*f_s;

subplot(2,2,3);
plot(t,y_c);
xlim([0 4]);
xlabel('time t');
ylabel('y_C(t)');
subplot(2,2,4);
plot(t,sC);
xlim([0 4]);
xlabel('time t');
ylabel('s_C(t)');


