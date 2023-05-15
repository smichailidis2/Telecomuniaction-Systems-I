function X = bits_to_2PAM(b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = bits_to_2PAM(b)                                                           %
%                                                                               %
% OUTPUT                                                                        %    
%       X : sequence of mapped symbols using 2-PAM mapping as shown below:      %
%                                                                               %
%           0 ---> + 1                                                          %
%           1 ---> - 1                                                          %
%                                                                               %
% INPUT                                                                         %      
%       b : input sequence of bits quantized to 0 & 1                           %
%                                                                               %
%    S. K. Michailidis, March 2023                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X = zeros(length(b));

% if b == 1, out <= -1*b elseif b == 0, out <= b + 1

X = -b .* ( b == 1 ) +  ( b + 1 ) .* ( b == 0 );
