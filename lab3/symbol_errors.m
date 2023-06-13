function num_of_symbol_errors = symbol_errors(estimated_X,actual_X)

Le = length(estimated_X);
num_of_symbol_errors = 0;

for i = 1 : Le
    if ( estimated_X(i,1) ~= actual_X(i,1) || estimated_X(i,2) ~= actual_X(i,2))
        num_of_symbol_errors = num_of_symbol_errors + 1;
    end
end


end