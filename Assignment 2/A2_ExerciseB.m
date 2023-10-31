close all
clear all
clc

    %Detection of a pulse: Exercise 2.b
    
    %defining variables

f_s = 100;               %sampling frequency
T_s = 1/f_s;             %sampling time
T = 1;                   %window time
ta = 0 : T_s : T-T_s;    %time axis
N = length(ta);          %total number of samples

f1 = 1;                 %frequency value 
A1 = sqrt(2);           %amplitude value 

s = A1*sin(2*pi*f1*ta);       %generating signal

E_s = sum(s.^2)*T_s;          %energy of the signal
P_s = (A1^2)/2;               %power of the signal using the theoretical formula

                              %energy and power are equal to 1
SNR = -10;               %SNR in dB
SNR_n = 10^(SNR/10);     %SNR in natural units = 0.100

P_n = P_s/SNR_n;         %noise power

n = sqrt(P_n) * randn(1,N);     %noise signal realization

r = s + n;               %generating signal plus noise


    %plotting 
    
figure
plot(ta,r);
title(sprintf('Example of signal plus noise\nSNR = %.2f dB', SNR));
xlabel('time');
ylabel('r');
grid on
hold on
plot(ta,s);
legend('r(t)','s(t)');


    %generate a number of noise realizations
    
num_sim = 100000;                    %number or simulations

n_r = sqrt(P_n) * randn(num_sim,N);  %num_sim noise signals
                                     
n_norm = zeros(num_sim,N);
r_norm = zeros(num_sim,N);


r_r(1:100) = n_r(1:100)+s;            %first noise plus signal
n_norm(1:100) = n_r(1:N)/(sqrt(sum(n_r(1:N).^2)/N));    %first noise signal normalized
r_norm(1:100) = r_r(1:N)/(sqrt(sum(r_r(1:N).^2)/N));    %first noise+signal normalized

for i = 1:num_sim-1   %normalize all other r and n
    
    tmp = i*N+1;
    r_r(tmp:tmp+N-1) = n_r(tmp:tmp+N-1)+s;
    n_norm(tmp:tmp+N-1) = n_r(tmp:tmp+N-1)/(sqrt(sum(n_r(tmp:tmp+N-1).^2)/N));
    r_norm(tmp:tmp+N-1) = r_r(tmp:tmp+N-1)/(sqrt(sum(r_r(tmp:tmp+N-1).^2)/N));
    
    %en = sum(n_norm(tmp:tmp+N-1).^2)*T_s
    %en = sum(r_norm(tmp:tmp+N-1).^2)*T_s;  <--  %uncomment to prove r and n
                                                 %are normalized
end

H0 = zeros(num_sim,1);
H1 = zeros(num_sim,1);

H0(1) = abs(sum(n_norm(1:N).*s)*T_s);       %correlation coeffs for first (n',s)
H1(1) = abs(sum(r_norm(1:N).*s)*T_s);       %correlation coeffs for first (r',s)

for x = 2:num_sim     %correlation coefficients for all other n' and r'
    
    temp = (x-1)*N+1;
    H0(x) = abs(sum(n_norm(temp:temp+N-1).*s)*T_s);  %correlation coeffs for (n',s)
    H1(x) = abs(sum(r_norm(temp:temp+N-1).*s)*T_s);  %correlation coeffs for (r',s)
    
end   


    %plotting

axis = 1:1:num_sim;
figure
plot(axis, H0)
grid on
hold on
plot(axis, H1)
title(sprintf('Correlation: simulation outcomes\nSNR = %.2f dB', SNR));
ylim([0 1]);
xlim([0 num_sim]);
xlabel('number of simulation');
ylabel('Γ');
legend('H0','H1');


    %histograms
   
figure;
histogram(H0,'Normalization','Probability')
hold on
grid on
histogram(H1,'Normalization','Probability')
title(sprintf('Correlation: pdf of Γ under H_0 and H_1\nSNR = %.2f dB', SNR));
legend('H0','H1');


    %threshold values
    
tmin = min(H0);
tmax = max(H1);
N_bins = 100;
t = linspace(tmax,tmin,N_bins);

P_fa = zeros(1,N_bins);
P_md = zeros(1,N_bins);

for z = 1:N_bins
    
    P_fa(z) = length(find(H0>=t(z)))/num_sim;       %false alarm probability
    P_md(z) = length(find(H1<t(z)))/num_sim;        %missed detection probability
    
    if(length(find(H1<t(z)))<30)        %for better reliability
        P_md(z) = 0;
    end
    
    if(length(find(H0>t(z)))<30)        %for better reliability
        P_fa(z) = 0;
    end
    
end


    %plotting
    
figure
semilogy(t, P_fa)
hold on
semilogy(t, P_md)
grid on
title(sprintf('Correlation: P_f_a and P_m_d vs. threshold\nSNR = %.2f dB', SNR));
legend('P_f_a','P_m_d');
xlabel('threshold t');
ylabel('P_f_a, P_m_d');


    %ROC curve
    
P_d = (1-P_md);

figure;
semilogx(P_fa,P_d);
loglog(P_fa,P_d);
title(sprintf('ROC curve\nSNR = %.2f dB', SNR));
grid on
xlabel('P_f_a');
ylabel('P_d');
ylim([0.29 1]);
legend('correlation ROC');


