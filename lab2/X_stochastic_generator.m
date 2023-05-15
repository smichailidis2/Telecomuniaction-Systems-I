function [X,t,T_total] = X_stochastic_generator(X_n,N,SRRC_pulse,srrc_time,T,over)
%X_STOHASTIC_GENERATOR Summary
%   Inputs :
%      X_n      --> symbols
%       N       --> number of symbols total in X_n
%      SRRC     --> truncated srrc pulse
%      over     --> oversampling factor of the SRC pulse
%   srrc_time   --> time vector of srrc pulse
%       T       --> symbol period
%
%   Output:
%       X    ---> Stochastic process
%       t    ---> time axis
%    T_total ---> total time duration of the Stochastic process

Ts = T/over;

x_delta = (1/Ts) * upsample(X_n, over);
x_delta_time=0:Ts:(N)*T-Ts;

X = conv(x_delta,SRRC_pulse)*Ts;

t = min(srrc_time) + min(x_delta_time) : Ts : max(srrc_time) + max(x_delta_time);
T_total = length(t)*Ts;

end

