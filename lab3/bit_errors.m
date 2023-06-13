function num_of_bit_errors = bit_errors(estimated_bit_sequence,actual_bit_sequence)

L = length(estimated_bit_sequence);
num_of_bit_errors = 0;

for i = 1 : L
    if  estimated_bit_sequence(i) ~= actual_bit_sequence(i)    
        num_of_bit_errors = num_of_bit_errors + 1;
    end
end


end