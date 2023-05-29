clear;
close all;
clc

%% Telecommunication Systems 1 %%
% EXERCISE 3
% Michailidis Stergios 
% 2020030080

%% PART 1 %%
%% 1)

N = 100;
b = (sign(randn(4*N, 1)) + 1)/2;

figure('Name','Question -1-','NumberTitle','off');
stem(b)
grid on;
axis([-1 401 -0.2 1.2])
title('Bit Train')

%% 2)
[Xi_n,Xq_n] = bits_to_16PSK(b);

figure('Name','Question -2-','NumberTitle','off');
stem(Xi_n,'r')
grid on
title('Mapped symbols (Real part)')
figure('Name','Question -2-','NumberTitle','off');
stem(Xq_n,'k')
grid on
title('Mapped symbols (Imaginary part)')

%% 3)

T = 10^(-2);
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 2;
a = 0.33;

[phi_t,t] = srrc_pulse(T,over,A,a);

% upsample the symbol trains
Xi_n_delta = (1/Ts) * upsample(Xi_n, over);
Xq_n_delta = (1/Ts) * upsample(Xq_n, over);
x_delta_time=0:Ts:(N)*T-Ts;

% cumpute the convoltuions
X_in = conv(phi_t,Xi_n_delta)*Ts;
X_Qn = conv(phi_t,Xq_n_delta)*Ts;

% new time axis
t_conv = min(t) + min(x_delta_time) : Ts : max(t) + max(x_delta_time);

% plots
figure('Name','Question -3-','NumberTitle','off');
plot(t_conv,X_in,'r')
grid on
xlabel('time')
ylabel('Amplitude of X_i_,_n ⋆ X_i_,_d_e_l_t_a')
title('Output Waveform')
figure('Name','Question -3-','NumberTitle','off');
plot(t_conv,X_Qn,'k')
grid on
xlabel('time')
ylabel('Amplitude of X_Q_,_n ⋆ X_Q_,_d_e_l_t_a')
title('Output Waveform')

% periodograms:
T_total = max(t_conv) - min(t_conv);
Nf = 2048;
F_axis = -Fs/2 : Fs/Nf : Fs/2 - Fs/Nf;

X_in_freq = fftshift(fft(X_in,Nf)*Ts);
X_Qn_freq = fftshift(fft(X_Qn,Nf)*Ts);

P_X_in = ( abs(X_in_freq).^2 ) / T_total;
P_X_Qn = ( abs(X_Qn_freq).^2 ) / T_total;

% plot the periodograms
figure('Name','Question -3-','NumberTitle','off');
plot(F_axis,P_X_in,'r')
grid on
xlabel('frequency')
ylabel('Periodogram of X_i_,_n ⋆ X_i_,_d_e_l_t_a')
title('Output Periodogram')
figure('Name','Question -3-','NumberTitle','off');
plot(F_axis,P_X_Qn,'k')
grid on
xlabel('frequency')
ylabel('Periodogram of X_Q_,_n ⋆ X_Q_,_d_e_l_t_a')
title('Output Periodogram')


%% 4)

F0 = 200;
cosinus_carrier =  2*cos(2*pi*F0*t_conv);
sinus_carrier   = -2*sin(2*pi*F0*t_conv);

X_I_t = X_in*cosinus_carrier;
X_Q_t = X_Qn*sinus_carrier;

% plots
figure('Name','Question -4-','NumberTitle','off');
plot(t_conv,X_I_t,'r')
grid on
xlabel('time')
ylabel('Amplitude of X_I')
title('Output Waveform')
figure('Name','Question -4-','NumberTitle','off');
plot(t_conv,X_Q_t,'k')
grid on
xlabel('time')
ylabel('Amplitude of X_Q')
title('Output Waveform')

% frequency domain
X_I_F = fftshift(fft(X_I_t,Nf)*Ts);
X_Q_F = fftshift(fft(X_Q_t,Nf)*Ts);

P_X_I = ( abs(X_I_F).^2 ) / T_total;
P_X_Q = ( abs(X_Q_F).^2 ) / T_total;

% plot the periodograms
figure('Name','Question -4-','NumberTitle','off');
plot(F_axis,P_X_I,'r')
grid on
xlabel('frequency')
ylabel('Periodogram of X_I')
title('Output Periodogram')
figure('Name','Question -4-','NumberTitle','off');
plot(F_axis,P_X_Q,'k')
grid on
xlabel('frequency')
ylabel('Periodogram of X_Q')
title('Output Periodogram')
