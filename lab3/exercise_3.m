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
X = bits_to_16PSK(b);

% euler's formula
Xi_n = real(X);
Xq_n = imag(X);


figure('Name','Question -2-','NumberTitle','off');
stem(Xi_n,'r')
axis([0 100 -1.2 1.2]);
grid on
title('Mapped symbols (Real part)')
figure('Name','Question -2-','NumberTitle','off');
stem(Xq_n,'k')
axis([0 100 -1.2 1.2]);
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
% The channel is ideal, h(t) = d_dirac(t)

%% 7)

% adding wgn
SNR_db = 20;

var_W = 1/ ( Ts * 10^(SNR_db/10) );

W_t = sqrt(var_W)*randn(1,length(t_conv))';

% channel output signal
Y_t = X_t + W_t;

%% 8)

% multiply the channel output

Y_i_t = Y_t.*cosinus_carrier*1/2;

Y_Q_t = Y_t.*sinus_carrier*1/2;

% time domain plots
figure('Name','Question -8.1-','NumberTitle','off');
plot(t_conv,Y_i_t,'r')
grid on
xlabel('time')
ylabel('Amplitude of Y_i')
title('Output Waveform')
figure('Name','Question -8.1-','NumberTitle','off');
plot(t_conv,Y_Q_t,'k')
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

% Y will be a complex vector with dim: 100x1. To make it compatible with
% the function detect_16_PSK.m , transform it to have dim: 100x2
Y = zeros(100,2);


% sampling the outputs

% We want to end up wtih N/4 = 100 samples at the output. To achieve this, start
% sampling from the samples 0 to 1000 with a sampling period equal to over = 10.

i = 1;
for j = 2*A*over + 1 : over : length(t_conv_new) - 2*A*over
   Y(i,1)=Y_i_filtered_t(j);
   Y(i,2)=Y_Q_filtered_t(j);
   i=i+1;
   assert(i<=101)
end

scatfig = scatterplot(Y);
scatfig.Name = 'Question -10-';
scatfig.NumberTitle = 'off';
title('Output Sequence Y');
grid on;


%% 11)

% detect the symbols from Y and estimate input bit sequence
[est_X,est_bit_seq] = detect_16_PSK(Y);


%% 12)

% calculate symbol errors
X_actual(:,1) = Xi_n;
X_actual(:,2) = Xq_n;

num_of_symbol_errors = symbol_errors(est_X,X_actual);

disp('Number of symbol errors:')
disp(num_of_symbol_errors)

%% 13)

% calculate bit errors
num_of_bit_errors = bit_errors(est_bit_seq',b);

disp('Number of bit errors:')
disp(num_of_bit_errors)


%% PART 2 %%
%% 1)

r = 1;

P_err_symbol = zeros(1,14);
P_err_bit    = zeros(1,14);

p_error_symbol_upper_bound = zeros(1,14);
p_error_bit_lower_bound    = zeros(1,14);

for SNRdb = -2:2:24
        
    total_number_of_symbol_errors = 0;
    total_number_of_bit_errors    = 0;


    for K = 1 : 1000
        % input bit seq
        b = (sign(randn(4*N, 1)) + 1)/2;
        % map to 16-PSK
        X = bits_to_16PSK(b);
        
        Xi_n = real(X);
        Xq_n = imag(X);
        
        % initialize srrc 
        [phi_t,t] = srrc_pulse(T,over,A,a);

        % upsample the symbol trains
        Xi_n_delta = (1/Ts) * upsample(Xi_n, over);
        Xq_n_delta = (1/Ts) * upsample(Xq_n, over);
        x_delta_time=0:Ts:(N)*T-Ts;

        % cumpute the convoltuions
        X_in = conv(phi_t,Xi_n_delta)*Ts;
        X_Qn = conv(phi_t,Xq_n_delta)*Ts;
        t_conv = min(t) + min(x_delta_time) : Ts : max(t) + max(x_delta_time);

        % multiply with carriers
        cosinus_carrier =  2*cos(2*pi*F0*t_conv)';
        sinus_carrier   = -2*sin(2*pi*F0*t_conv)';
        
        X_I_t = X_in.*cosinus_carrier;
        X_Q_t = X_Qn.*sinus_carrier;

        % transmitter output 
        X_t = X_I_t + X_Q_t;
        var_W = 1/ ( Ts * 10^(SNRdb/10) );
        W_t = sqrt(var_W)*randn(1,length(t_conv))';

        % channel output signal
        Y_t = X_t + W_t;

        % multiply with carriers
        Y_i_t = Y_t.*cosinus_carrier*1/2;
        Y_Q_t = Y_t.*sinus_carrier*1/2;

        % cumpute the final convoltuions
        Y_i_filtered_t = conv(Y_i_t,phi_t)*Ts;
        Y_Q_filtered_t = conv(Y_Q_t,phi_t)*Ts;

        % sampling
        Y = zeros(100,2);
        u = 1;
        for p = 2*A*over + 1 : over : length(t_conv_new) - 2*A*over
           Y(u,1)=Y_i_filtered_t(p);
           Y(u,2)=Y_Q_filtered_t(p);
           u=u+1;
           assert(u<=101)
        end
        

        % detect the symbols from Y and estimate input bit sequence
        [est_X,est_bit_seq] = detect_16_PSK(Y);

        % Total symbol errors
        X_actual(:,1) = Xi_n;
        X_actual(:,2) = Xq_n;
        
        se = symbol_errors(est_X,X_actual);
        total_number_of_symbol_errors = total_number_of_symbol_errors + se;

        % Total bit errors
        be = bit_errors(est_bit_seq',b);
        total_number_of_bit_errors = total_number_of_bit_errors + be;

    end
% experimental symbol and bit errors
P_err_symbol(r) = (total_number_of_symbol_errors)/(K*length(X));
P_err_bit(r)    = (total_number_of_bit_errors)/(K*length(b));


% translate SNRdb to SNR : SNRdb = 10*log_10(SNR) => SNR = 10^(SNRdb/10)
SNR = 10^(SNRdb/10);

% theoretical symbol,bit errors
p_error_symbol_upper_bound(r) = 2 * Q( sqrt( 2*(SNR) )*sin( pi/16 ) );
p_error_bit_lower_bound(r) = p_error_symbol_upper_bound(r)/4;

r = r + 1;

end

SNRdb_axis = -2 : 2 : 24;

figure('Name','Question -B.2-','NumberTitle','off');
semilogy(SNRdb_axis,p_error_symbol_upper_bound,'k','LineWidth',2)
hold on;
semilogy(SNRdb_axis,P_err_symbol,'r')
grid on;
xlabel('SNR_d_B')
ylabel('Symbol error Probability')
legend('Upper Bound','Total Symbol error')
title('Total Symbol error as function of SNR_d_B')
hold off


figure('Name','Question -B.3-','NumberTitle','off');
semilogy(SNRdb_axis,p_error_bit_lower_bound,'k','LineWidth',2)
hold on;
semilogy(SNRdb_axis,P_err_bit,'r')
grid on;
xlabel('SNR_d_B')
ylabel('Bit error Probability')
legend('Lower Bound','Total Bit error')
title('Bit error as function of SNR_d_B')
hold off

fprintf('\n---------Part B---------\n')
fprintf('\nTotal Symbol error probability per SNR_d_B:\n\n')
disp('    SNRdb  | P_err_symbol')
disp([SNRdb_axis' P_err_symbol'])

fprintf('Total bit error probability per SNR_d_B:\n\n')
disp('    SNRdb  | P_err_bit')
disp([SNRdb_axis' P_err_bit'])

