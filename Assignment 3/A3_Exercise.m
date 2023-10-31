clear all
close all
clc


    %defining variables

f_s = 100;          %sampling frequency
T_s = 1/f_s;        %sampling time
T = 30;             %observation window
t = 0: T_s : T-T_s;         %time axis
T_skip = 10;                %time to skip
sample_start = T_skip*f_s+1;    %sample to start from
N = length(t);                  %number of samples to observe           

f_res = f_s/N;                  %frequency resolution
f = 0:f_res:f_s-f_res;          %frequency axis
f_symm = f-f_s/2;               %symmetric frequency axis


    %uploading the text file
    
filename = 'pulse.txt';
delimiterln = ' ';
headerlinesln = 1;
Data_struct = importdata(filename, delimiterln, headerlinesln);
Led_R_in = Data_struct.data(:,1);
Led_IR_in = Data_struct.data(:,2);

sample_size = length(Led_R_in);     %number of samples

Led_R = Led_R_in(sample_start:sample_size);     %resize to skip first 10 seconds
Led_IR = Led_IR_in(sample_start:sample_size);   %resize to skip first 10 seconds

Led_R_plot = Led_R(1:N);
Led_IR_plot = Led_IR(1:N);

F_Led_R = fft(Led_R_plot)*T_s;      %Fourier transform of the red signal
MR = abs(F_Led_R);
MR_shift = fftshift(MR);

F_Led_IR = fft(Led_IR_plot)*T_s;    %Fourier transform of the infrared signal
MIR = abs(F_Led_IR);
MIR_shift = fftshift(MIR);


    %low-pass filter

H_L = rectangularPulse(-3,3,f_symm);

Y_R = fftshift(F_Led_R).*transpose(H_L);   %low pass filter on R Fourier signal
Y_R_shift = ifftshift(Y_R);
y_r = ifft(Y_R_shift)*f_s;

Y_IR = fftshift(F_Led_IR).*transpose(H_L); %low pass filter on IR Fourier signal
Y_IR_shift = ifftshift(Y_IR);
y_ir = ifft(Y_IR_shift)*f_s;

    %plotting
    
figure(1)
subplot(2,2,1);
plot(t,Led_R_plot,'LineWidth',0.1);
grid on
xlabel('Time [s]');
ylabel('Led_R');
title(sprintf('Original Red signal\n'));
subplot(2,2,2);
plot(t,Led_IR_plot,'LineWidth',0.1);
grid on
xlabel('Time [s]');
ylabel('Led_I_R');
title(sprintf('Original Infrared signal\n'));

subplot(2,2,3);
plot(t,y_r,'LineWidth',0.1);
grid on
xlabel('Time [s]')
ylabel('Led_R')
title(sprintf('Red signal after low pass\n'));
subplot(2,2,4);
plot(t,y_ir,'LineWidth',0.1);
grid on
xlabel('Time [s]')
ylabel('Led_I_R')
title(sprintf('Infrared signal after low pass\n'));


    %finding min and max and applying interpolation

[pks1_r,locs1_r]=findpeaks(y_r,t);
[pks2_r,locs2_r]=findpeaks(-y_r,t);
pks2_r=-pks2_r;
HH1_r=interp1(locs1_r,pks1_r,t,'spline');
HH2_r=interp1(locs2_r,pks2_r,t,'spline');

[pks1_ir,locs1_ir]=findpeaks(y_ir,t);
[pks2_ir,locs2_ir]=findpeaks(-y_ir,t);
pks2_ir=-pks2_ir;
HH1_ir=interp1(locs1_ir,pks1_ir,t,'spline');
HH2_ir=interp1(locs2_ir,pks2_ir,t,'spline');

    
    %computing R
    
R = ((HH1_r-HH2_r)./HH2_r)./((HH1_ir-HH2_ir)./HH2_ir);
R_mean = mean(R);
SaO2 = 110 - 25*R_mean;

figure
subplot(2,1,1)
plot(t,y_r);
hold on
plot(t,HH1_r,'r')
hold on
plot(t,HH2_r,'r')
grid on
title(sprintf('RED'));

subplot(2,1,2)
plot(t,y_ir);
hold on
plot(t,HH1_ir,'r')
hold on
plot(t,HH2_ir,'r')
grid on
title(sprintf('INFRARED SPO_2 = %.2f%%', SaO2));
   

    %high-pass filter
 
H_H1 = rectangularPulse(-0.5,0.5,f_symm); 
H_H = 1-H_H1;
    
F_Led_R2 = fft(y_r)*T_s;
Y_R2 = fftshift(F_Led_R2).*transpose(H_H);
Y_R2_shift = ifftshift(Y_R2);
y_r2 = ifft(Y_R2_shift)*f_s;  %high and low pass filtered RED signal


    %pulse rate computation
    
[pks1_r2,locs1_r2]=findpeaks(y_r2,t);
max_d = diff(locs1_r2);
d_mean = mean(max_d);
BPM = 60/d_mean;

figure(3)
findpeaks(y_r2,t)
grid on
xlabel('Time [s]');
xlim([0 30]);
title(sprintf('filtered red signal: BPM = %.2f', BPM));


    %alternative pulse rate computation
    
F_RED = fft(y_r2)*T_s;
M_RED = abs(F_RED);
M_shifted = fftshift(M_RED);

BPM_axis = f_symm*60;
[p,l] = findpeaks(M_shifted,BPM_axis);
len = length(p);
p = p(len/2:len);
l = l(len/2:len);

[p_max, index] = max(p);
l_max = l(index);

figure
plot(BPM_axis,M_shifted)
grid on
xlabel('BPM')
xlim([60 90])

title(sprintf('filtered red signal - frequency axis - BPM = %.2f',l_max));

