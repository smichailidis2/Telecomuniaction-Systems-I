function [est_X,est_bit_seq] = detect_16_PSK(Y)

m = 16; % 16-psk

le = length(Y);
est_X = zeros(2,le)';
est_bit_seq = zeros(1,4*le);

% nearest neighbour algorithm ---------------------------------------------
distances = zeros(1,16);

for i = 1:le

    for v = 1:m

        % Calculate distance from each symbol
        const_symb = exp(2*pi*1i*(v-1)*1/m);

        distances(v) = sqrt( ( Y(i,1) - real(const_symb) ).^2 + ( Y(i,2) - imag(const_symb) ).^2 );
    end
    

    % Nearest neighbour for the i-th symbol
    nearest_nb = min(distances);
    
    % Estimation 
    for n = 1:m
        estim_symb = exp(2*pi*1i*(n-1)*1/m);
        %disp(estim_symb)

        if( distances(n) == nearest_nb )
            est_X(i,1) = real(estim_symb);
            est_X(i,2) = imag(estim_symb);
        end

    end
   
end

%disp(size(est_X))

% -------------------------------------------------------------------------

% Execute the inverse procedure from the "bits_to_16PSK.m"

c = 1;

M = 0:15;

symbols = exp(2*pi*1i*M*1/16);
%disp(symbols)

for k = 1 : 4 : 4*le-3
    
    % est_X = 1 ---> '0000'
    if(est_X(c,1) == real(symbols(1)) && est_X(c,2) == imag(symbols(1)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 0;
    end
    

    % est_X = exp(2*pi*1i*1/16) ---> '1000'
    if(est_X(c,1) == real(symbols(2)) && est_X(c,2) == imag(symbols(2)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*2/16) ---> '1001'
    if(est_X(c,1) == real(symbols(3)) && est_X(c,2) == imag(symbols(3)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 1;
    end


    % est_X = exp(2*pi*1i*3/16) ---> '1011'
    if(est_X(c,1) == real(symbols(4)) && est_X(c,2) == imag(symbols(4)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 1;
    end


    % est_X = exp(2*pi*1i*4/16) ---> '1010'
    if(est_X(c,1) == real(symbols(5)) && est_X(c,2) == imag(symbols(5)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*5/16) ---> '1110'
    if(est_X(c,1) == real(symbols(6)) && est_X(c,2) == imag(symbols(6)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*6/16) ---> '1111'
    if(est_X(c,1) == real(symbols(7)) && est_X(c,2) == imag(symbols(7)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 1;
    end


    % est_X = exp(2*pi*1i*7/16) ---> '1101'
    if(est_X(c,1) == real(symbols(8)) && est_X(c,2) == imag(symbols(8)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 1;
    end
    

    % est_X = exp(2*pi*1i*8/16) ---> '1100'
    if(est_X(c,1) == real(symbols(9)) && est_X(c,2) == imag(symbols(9)))
        est_bit_seq(k)   = 1;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*9/16) ---> '0100'
    if(est_X(c,1) == real(symbols(10)) && est_X(c,2) == imag(symbols(10)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*10/16) ---> '0101'
    if(est_X(c,1) == real(symbols(11)) && est_X(c,2) == imag(symbols(11)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 1;
    end

    % est_X = exp(2*pi*1i*11/16) ---> '0111'
    if(est_X(c,1) == real(symbols(12)) && est_X(c,2) == imag(symbols(12)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 1;
    end

    % est_X = exp(2*pi*1i*12/16) ---> '0110'
    if(est_X(c,1) == real(symbols(13)) && est_X(c,2) == imag(symbols(13)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 1;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*13/16) ---> '0010'
    if(est_X(c,1) == real(symbols(14)) && est_X(c,2) == imag(symbols(14)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 0;
    end


    % est_X = exp(2*pi*1i*14/16) ---> '0011'
    if(est_X(c,1) == real(symbols(15)) && est_X(c,2) == imag(symbols(15)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 1;
        est_bit_seq(k+3) = 1;
    end


    % est_X = exp(2*pi*1i*15/16) ---> '0001'
    if(est_X(c,1) == real(symbols(16)) && est_X(c,2) == imag(symbols(16)))
        est_bit_seq(k)   = 0;
        est_bit_seq(k+1) = 0;
        est_bit_seq(k+2) = 0;
        est_bit_seq(k+3) = 1;
    end
    %disp(c)
    c = c + 1;
    assert(c<=101)
end


end
