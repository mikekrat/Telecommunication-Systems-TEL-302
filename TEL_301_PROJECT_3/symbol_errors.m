function num_of_symbol_errors=symbol_errors(est_X,X)
num_of_symbol_errors=0;
for i=1:length(X)
    if(abs(X(i,1)-est_X(i,1))>0.001)
        num_of_symbol_errors=num_of_symbol_errors+1;
    end
end
for i=1:length(X)
    if(abs(X(i,2)-est_X(i,2))>0.001)
        num_of_symbol_errors=num_of_symbol_errors+1;
    end
end

num_of_symbol_errors=ceil(num_of_symbol_errors/2);

end