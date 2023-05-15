function Xn = bits_to_4PAM(b1,b2)

% X = bits_to_4PAM(b1,b2)                                                       
%                                                                               
% OUTPUT                                                                        
%       Xn : sequence of mapped symbols using 4-PAM mapping as shown below:     
%                                                                              
%           00 ---> + 3                                                         
%           01 ---> + 1                                                     
%           11 ---> - 1
%           10 ---> - 3                                                         
%                                                                               
% INPUT                                                                         
%       b1,2 : input sequences of bits quantized to 0 & 1                           
%                                                                               
%    S. K. Michailidis, May 2023                                              

% Xn = zeros(length(b1));

Xn = -3 * b1 .* (b1 == 1 & b2 == 0) - 1 * b1 .* (b1 == 1 & b2 == 1) + b2  .* (b1 == 0 & b2 == 1) + 3 * (b2 + 1) .* (b1 == 0 & b2 ==0);

