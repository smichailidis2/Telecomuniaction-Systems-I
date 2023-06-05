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
cosinus_carrier =  2*cos(2*pi*F0*t_conv)';
sinus_carrier   = -2*sin(2*pi*F0*t_conv)';

X_I_t = X_in.*cosinus_carrier;
X_Q_t = X_Qn.*sinus_carrier;

% plots
figure('Name','Question -4.1-','NumberTitle','off');
plot(t_conv,X_I_t,'r')
grid on
xlabel('time')
ylabel('Amplitude of X_I')
title('Output Waveform')
figure('Name','Question -4.1-','NumberTitle','off');
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
figure('Name','Question -4.2-','NumberTitle','off');
plot(F_axis,P_X_I,'r')
grid on
xlabel('frequency')
ylabel('Periodogram of X_I')
title('Output Periodogram')
figure('Name','Question -4.2-','NumberTitle','off');
plot(F_axis,P_X_Q,'k')
grid on
xlabel('frequency')
ylabel('Periodogram of X_Q')
title('Output Periodogram')


%% 5)

% transmitter output - time domain
X_t = X_I_t + X_Q_t;

figure('Name','Question -5-','NumberTitle','off');
plot(t_conv,X_t);
grid on;
xlabel('time')
ylabel('Amplitude of X')
title('Transmitter Output Waveform')

% transmitter output - frequency domain
X_F = X_I_F + X_Q_F;    % fft is a linear operator

P_X = abs(X_F).^2/T_total;

figure('Name','Question -5-','NumberTitle','off');
plot(F_axis,P_X)
grid on
xlabel('frequency')
ylabel('Periodogram of X')
title('Transmitter Output Periodogram')


%% 6)
% The channel is ideal

%% 7)

% adding wgn
SNR_db = 20;

var_W = 1/ ( Ts * 10^(SNR_db/10) );
var_N = Ts*var_W/2;

W_t = sqrt(var_N)*randn(1,length(t_conv))';

% channel output signal
Y_t = X_t + W_t;

%% 8)

% multiply the channel output

Y_i_t = Y_t.*cosinus_carrier;

Y_Q_t = Y_t.*sinus_carrier;

% time domain plots
figure('Name','Question -8.1-','NumberTitle','off');
plot(t_conv,X_in,'r')
grid on
xlabel('time')
ylabel('Amplitude of Y_i')
title('Output Waveform')
figure('Name','Question -8.1-','NumberTitle','off');
plot(t_conv,X_Qn,'k')
grid on
xlabel('time')
ylabel('Amplitude of Y_Q')
title('Output Waveform')


% frequency domain

Y_i_F = fftshift(fft(Y_i_t,Nf)*Ts);

Y_Q_F = fftshift(fft(Y_Q_t,Nf)*Ts);


P_Y_i = ( abs(Y_i_F).^2 ) / T_total;
P_Y_Q = ( abs(Y_Q_F).^2 ) / T_total;

% plot the periodograms
figure('Name','Question -8.2-','NumberTitle','off');
plot(F_axis,P_Y_i,'r')
grid on
xlabel('frequency')
ylabel('Periodogram of Y_i')
title('Output Periodogram')
figure('Name','Question -8.2-','NumberTitle','off');
plot(F_axis,P_Y_Q,'k')
grid on
xlabel('frequency')
ylabel('Periodogram of Y_Q')
title('Output Periodogram')


%% 9)

Y_i_filtered_t = conv(Y_i_t,phi_t)*Ts;

Y_Q_filtered_t = conv(Y_Q_t,phi_t)*Ts;

% calculate new time vector and new T_total
t_conv_new = min(t_conv) + min(t) : Ts : max(t_conv) + max(t);
T_total_new = max(t_conv_new) - min(t_conv_new);


% time domain plots
figure('Name','Question -9.1-','NumberTitle','off');
plot(t_conv_new,Y_i_filtered_t,'r')
grid on
xlabel('time')
ylabel('Amplitude of Y_i_ _f_i_l_t_e_r_e_d')
title('Output Waveform')
figure('Name','Question -9.1-','NumberTitle','off');
plot(t_conv_new,Y_Q_filtered_t,'k')
grid on
xlabel('time')
ylabel('Amplitude of Y_Q_ _f_i_l_t_e_r_e_d')
title('Output Waveform')


% frequency domain

Y_i_filtered_F = fftshift(fft(Y_i_filtered_t,Nf)*Ts);

Y_Q_filtered_F = fftshift(fft(Y_Q_filtered_t,Nf)*Ts);


P_Y_i_filtered = ( abs(Y_i_filtered_F).^2 ) / T_total_new;
P_Y_Q_filtered = ( abs(Y_Q_filtered_F).^2 ) / T_total_new;

% plot the periodograms
figure('Name','Question -9.2-','NumberTitle','off');
plot(F_axis,P_Y_i_filtered,'r')
grid on
xlabel('frequency')
ylabel('Periodogram of Y_i_ _f_i_l_t_e_r_e_d')
title('Output Periodogram')
figure('Name','Question -9.2-','NumberTitle','off');
plot(F_axis,P_Y_Q_filtered,'k')
grid on
xlabel('frequency')
ylabel('Periodogram of Y_Q_ _f_i_l_t_e_r_e_d')
title('Output Periodogram')


%% 10)

% TODO (to explain)

Y = zeros(100,2);

i = 1;

for j = 2*A*over + 1 : over : length(t_conv_new) - 2*A*over
   Y(i,1)=Y_i_filtered_t(j);
   Y(i,2)=Y_Q_filtered_t(j);
   i=i+1;
end

scatterplot(Y);
grid on;
title('samples');

