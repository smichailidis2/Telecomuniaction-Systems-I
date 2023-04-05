
% ======= FIRST EXERCISE THL 1 ======= %
% Michailidis Stergios, 2020030080
% March 2023

% ---Question A---

clear;
close all;
clc

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
% time vector t is the same for all phis...

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

f_axis = -Fs/2 : Fs/N_f : Fs/2 - Fs/N_f;

% Fourrier transforms
PHI_1 = fftshift( fft(phi1,N_f) )*Ts;
PHI_2 = fftshift( fft(phi2,N_f) )*Ts;
PHI_3 = fftshift( fft(phi3,N_f) )*Ts;

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

% Constant function, hfactor can be tweaked to approximate most optimal
% bandwidth
hfactor = 10^3;
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
cstr = sprintf('C = %g ', T/hfactor);

grid on;
title('Question A-3')
xlabel('Frequency Spectrum')
ylabel('P.S.D. of SRRC (Logarithmic Scale)')
legend([s1,s2,s3,s4],'a = 0','a = 0.5','a = 1',cstr)


% ---Question B---

clear;
close all;
clc

% ( B-1 )

% questions 1. & 2.

T = 10^(-3);
over = 10;
A = 4;
Ts = T/over;
Fs = 1/Ts;

for a = 0:0.5:1
    for k = 0:2
        [phi,t] = srrc_pulse(T,over,A,a);
        shifted_phi = [zeros(1,k*over) phi(1:end-k*over)];
        figure()
        plot(t,phi)
        hold on;
        plot(t,shifted_phi)
        hold on;
        hold off;
        grid on;
        title('SRRC pulses')
        
        % Hadamard product of the two functions
        prod = phi .* shifted_phi;
        figure()
        plot(t, prod)
        title('Hadamard product of the srrc pulses')
        grid on;

    end
end



% question 3.

clear;
close all;
clc

T = 10^(-3);
over = 10;
A = 4;
Ts = T/over;

for a = 0:0.5:1
    for k = 0:3
        [phi,t] = srrc_pulse(T,over,A,a);
        shifted_phi = [zeros(1,k*over) phi(1:end-k*over)];
        
        % Hadamard product of the two functions
        prod = phi .* shifted_phi;
        
        % simulating the integral
        integ = sum( prod )*Ts
        
        % plotting the integral values
        integ_plot = integ.*(t.^0);
        txt = ['Rff = ',num2str(integ)];
        plot(t,integ_plot,'DisplayName',txt)
        hold on;
    end
end
grid on;
axis([-4*T 4*T -0.2 1.2])
legend show
hold off;



% ---Question C---

clear;
close all;
clc

T = 10^(-3);
over = 10;
a = 0.5;
A = 4;
Ts = T/over;
Fs = 1/Ts;
N = 100;

% ( C-1 )

b = (sign(randn(N, 1)) + 1)/2;
figure
stem(b)

% ( C-2 )

%  a) 2-PAM :
X_k = bits_to_2PAM(b);
figure
stem(X_k);

% b) Simulation of the sum:

X_delta = 1/Ts * upsample(X_k, over);
time = 0:Ts:N*T-Ts;
figure
stem(time,X_delta);
grid on;

% c) Simulation of the first convolution:

[phi,t] = srrc_pulse(T,over,A,a);

% compute the convolution time
t_conv = min(time) + min(t) : Ts : max(time) + max(t) ;

X = conv(X_delta,phi)*Ts;
figure
plot(t_conv,X)
grid on;

% d) Simulation of the second convolution:

% First flip the srrc pulse and the time vectors:
phi_flipped = fliplr(phi);
t_new = fliplr(t);

% Compute the new convolution time
t_conv_new = min(t_conv) + min(t_new) : Ts : max(t_conv) + max(t_new) ;

Z = conv(X,phi_flipped)*Ts;

figure
plot(t_conv_new,Z)
grid on;
hold on;
stem([0:N-1]*T,X_k);
grid on;
