clear;
close all;
clc

% ======= FIRST EXERCISE THL 1 ======= %
% Michailidis Stergios, 2020030080
% March 2023

% ---Question A---

% ( A-1 )

T = 10^(-3);
over = 10;
A = 4;
Ts = T/over;
Fs = 1/Ts;

a1 = 0;
a2 = 0.5;
a3 = 1;

% Get the sqaure root raised cosine pulses
[phi1,~] = srrc_pulse(T,over,A,a1);
[phi2,~] = srrc_pulse(T,over,A,a2);
[phi3,t] = srrc_pulse(T,over,A,a3);

figure(1)
pl1 = plot(t,phi1,'k');
hold on;
pl2 = plot(t,phi2,'r');
hold on;
pl3 = plot(t,phi3,'b');

grid on;
axis([-4*T 4*T -10 50])
title('Question A-1')
xlabel('Time')
ylabel('Magnitude of SRRC')
legend([pl1,pl2,pl3],'a = 0','a = 0.5','a = 1')


% ( A-2 )

N_f = 2048;

f_axis = [-Fs/2 : Fs/N_f : Fs/2 - Fs/N_f];

% Fourrier transforms
PHI_1 = fftshift( fft(phi1,N_f) );
PHI_2 = fftshift( fft(phi2,N_f) );
PHI_3 = fftshift( fft(phi3,N_f) );

% Energy Spectral Density (E.S.D.)
esd1 = abs(PHI_1).^2;
esd2 = abs(PHI_2).^2;
esd3 = abs(PHI_3).^2;

figure(2)
pl4 = plot(f_axis,esd1,'k');
hold on;
pl5 = plot(f_axis,esd2,'r');
hold on;
pl6 = plot(f_axis,esd3,'b');

grid on;
title('Question A-2-a')
xlabel('Frequency Spectrum')
ylabel('E.S.D. of SRRC')
legend([pl4,pl5,pl6],'a = 0','a = 0.5','a = 1')


% logarithmic scale y-axis
figure(3)
s1 = semilogy(f_axis,esd1,'k');
hold on;
s2 = semilogy(f_axis,esd2,'r');
hold on;
s3 = semilogy(f_axis,esd3,'b');

grid on;
title('Question A-2-b')
xlabel('Frequency Spectrum')
ylabel('E.S.D. of SRRC (Logarithmic Scale)')
legend([s1,s2,s3],'a = 0','a = 0.5','a = 1')


% ( A-3 )

% Theoretical bandwidth for each roll-off factor value
BW1 = (1 + a1)/( 2*T )
BW2 = (1 + a2)/( 2*T )
BW3 = (1 + a3)/( 2*T )

% Constant function
hfactor = 10^-3;
c = T/hfactor*(f_axis.^0); % vector of magnitude T/Cvalue and length f_axis

figure(4)
s1 = semilogy(f_axis,esd1,'k');
hold on;
s2 = semilogy(f_axis,esd2,'r');
hold on;
s3 = semilogy(f_axis,esd3,'b');
hold on;
s4 = semilogy(f_axis,c,'g');

% for the legend
gravstr = sprintf('c = %g ', T/Cvalue);

grid on;
title('Question A-3')
xlabel('Frequency Spectrum')
ylabel('P.S.D. of SRRC (Logarithmic Scale)')
legend([s1,s2,s3,s4],'a = 0','a = 0.5','a = 1',gravstr)
