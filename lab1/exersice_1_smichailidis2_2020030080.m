clear;
close all;
clc

% ======= FIRST EXERCISE THL 1 ======= %
% Michailidis Stergios, 2020030080
% March 2023

% ---Question A---

% A-1)

T = 10^(-3);
over = 10;
A = 4;
Ts = T/over;

a1 = 0;
a2 = 0.5;
a3 = 1;

[phi1,t1] = srrc_pulse(T,over,A,a1);
[phi2,t2] = srrc_pulse(T,over,A,a2);
[phi3,t3] = srrc_pulse(T,over,A,a3);

figure(1)
pl1 = plot(t1,phi1,'k');
hold on;
pl2 = plot(t2,phi2,'r');
hold on;
pl3 = plot(t3,phi3,'b');

title('Question A-1')
xlabel('Time')
ylabel('Magnitude of SRRC')

legend([pl1,pl2,pl3],'a = 0','a = 0.5','a = 1')

% A-2)



