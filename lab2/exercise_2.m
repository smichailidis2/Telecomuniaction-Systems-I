clear;
close all;
clc

%% Telecommunication Systems 1 %%
% Michailidis Stergios 
% 2020030080


%% PART A %%
%% 1)
clear;
close all;
clc

Nf = 2048;
%Nf = 4096;

T = 10^(-2);
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 4;
a = 0.5;

% phi_t -> srrc pulse 
[phi_t, time ]= srrc_pulse(T,over,A,a);

% compute fourrier transform of phi_t
PHI_F  = fftshift(fft(phi_t,Nf))*Ts;
F_axis = -(Fs/2) : (Fs/Nf) : (Fs/2) - (Fs/Nf);

figure
semilogy(F_axis,abs(PHI_F).^2,'k')
xlabel('Frequency')
ylabel('Energy Spectral Density of phi_t in log scale')
grid on;
title("A-1")


%% 2)
N = 100;

% generate random bit train
b = (sign(randn(N, 1)) + 1)/2;

% 2 - PAM map
X_n = bits_to_2PAM(b);

% generate the stohastic process X
[X,time_new,T_total] = X_stochastic_generator(X_n,N,phi_t,time,T,over);

figure
plot(time_new,X,'b')
title('Stochastic Process X')
xlabel("Time")
ylabel("Time evolution of X")
grid on;
title("A-2")

%Compute the theoretical power spectral density (P.S.D.)
S_X_theoretical = (var(X_n)/T)*abs(PHI_F.^2);
figure
semilogy(F_axis,S_X_theoretical,'b')
xlabel("Frequency")
ylabel("Power Spectral Density (P.S.D.) of X in log scale");
grid on;
title("A-2")


%% 3)

X_F = fftshift(fft(X,Nf))*Ts;
Px_F = ( abs(X_F).^2 )/T_total;
figure
plot(F_axis,Px_F,'m')
xlabel('Frequency')
ylabel('Periodogram of X')
grid on;
title("A-3.1")

figure
semilogy(F_axis,Px_F,'b')
xlabel('Frequency')
ylabel('Periodogram of X in log scale')
grid on;
title("A-3.1")

% k = 2;
% k = 10;
% k = 20;
 k = 50;
% k = 100;
% k = 500;
% k = 2000;

%  Let S_X_hat be the estimated P.S.D.
% We will generate k stochastic processes X_j.
% The adder matrix will have rows corresponding to the periodogram of each
% realization of each stochastic process generated. Afterwards, we shall
% take the mean of the matrix will result to the row vector S_X_hat.

Adder_matrix = zeros(k,Nf);


for j = 1:k
    % generate random bit train
    b = (sign(randn(N, 1)) + 1)/2;
    
    % 2 - PAM map
    X_n = bits_to_2PAM(b);

    [X_j,~,total_dur] = X_stochastic_generator(X_n,N,phi_t,time,T,over);

    %compute the periodograms
    X_Fj = fftshift(fft(X_j,Nf))*Ts;
    P_x_Fj = abs(X_Fj.^2)/total_dur;

    Adder_matrix(j,:) = P_x_Fj;
end

S_X_hat = mean(Adder_matrix);

figure;
semilogy(F_axis,S_X_hat,'r');
grid on;
hold on
semilogy(F_axis,S_X_theoretical,'b');
grid on;
legend('Experimental P.S.D.','Theoretical P.S.D.');
xlabel('Frequency')
ylabel('P.S.D. in log scale')
title("A-3")

%% 4)

% k = 2;
% k = 10;
% k = 20;
 k = 50;
% k = 100;
% k = 500;
% k = 2000;

Adder_matrix = zeros(k,Nf);


for j = 1:k
    % generate random bit train
    b1 = (sign(randn(N/2, 1)) + 1)/2;
    b2 = (sign(randn(N/2, 1)) + 1)/2;
    
    % 4 - PAM map
    X_n = bits_to_4PAM(b1,b2);

    [X_j,~,total_dur] = X_stochastic_generator(X_n,N/2,phi_t,time,T,over);

    %compute the periodograms
    X_Fj = fftshift(fft(X_j,Nf))*Ts;
    P_x_Fj = abs(X_Fj.^2)/total_dur;

    Adder_matrix(j,:) = P_x_Fj;
end

S_X_hat = mean(Adder_matrix);
S_X_theoretical = (var(X_n)/T)*abs(PHI_F.^2);


figure;
semilogy(F_axis,S_X_hat,'r');
grid on;
hold on
semilogy(F_axis,S_X_theoretical,'b');
grid on;
legend('Experimental P.S.D.','Theoretical P.S.D.');
xlabel('Frequency')
ylabel('P.S.D. in log scale')
title("A-4")

%% 5)

T_new = 2*T;
over_new = 2*over;
% Ts_new = Ts_old = Ts <=>  Fs_new = Fs_old = Fs

% srrc pulse and fourrier
[phi,t] = srrc_pulse(T_new,over_new,A,a);
PHI_F  = fftshift(fft(phi,Nf))*Ts;

N = 100;

% generate random bit train
b = (sign(randn(N, 1)) + 1)/2;

% 2 - PAM map
X_n = bits_to_2PAM(b);

% generate the stohastic process X
[X,st_time,T_total] = X_stochastic_generator(X_n,N,phi,t,T_new,over_new);

%Compute the theoretical power spectral density (P.S.D.)
S_X_theoretical = (var(X_n)/T_new)*abs(PHI_F.^2);

% Compute the fourrier transform of the stochastic process X and the
% periodogram
X_F = fftshift(fft(X,Nf))*Ts;
Px_F = ( abs(X_F).^2 )/T_total;
figure
plot(F_axis,Px_F,'m')
xlabel('Frequency')
ylabel('Periodogram of X')
grid on;
title("A-5")

figure
semilogy(F_axis,Px_F,'b')
xlabel('Frequency')
ylabel('Periodogram of X in log scale')
grid on;
title("A-5")

% same process as A_3 ...
% k = 2;
% k = 10;
% k = 20;
 k = 50;
% k = 100;
% k = 500;
% k = 2000;
Adder_matrix = zeros(k,Nf);
for j = 1:k
    % generate random bit train
    b = (sign(randn(N, 1)) + 1)/2;
    % 2 - PAM map
    X_n = bits_to_2PAM(b);

    [X_j,~,total_dur] = X_stochastic_generator(X_n,N,phi,t,T_new,over_new);
    %compute the periodograms
    X_Fj = fftshift(fft(X_j,Nf))*Ts;
    P_x_Fj = abs(X_Fj.^2)/total_dur;

    Adder_matrix(j,:) = P_x_Fj;
end

S_X_hat = mean(Adder_matrix);

figure;
semilogy(F_axis,S_X_hat,'r');
grid on;
hold on
semilogy(F_axis,S_X_theoretical,'b');
grid on;
legend('Experimental P.S.D.','Theoretical P.S.D.');
xlabel('Frequency')
ylabel('P.S.D. in log scale')
title("A-5")

%% PART B %%

F0 = 200;
t = linspace(-A,A,2*F0);

X = randn(1,5);         % normal Gaussian vector 
PHI = (2*pi)*rand(1,5); % uniform vector

figure
hold on;
for i = 1 : 5
    Y = X(i)*cos(2*pi*F0*t + PHI(i));
    txt = ['Y',num2str(i)];
    plot(t,Y,'DisplayName',txt)
end
grid on;
hold off
xlabel("time axis")
ylabel("Realizations of the stochastic processes")
legend show;
title("B")


