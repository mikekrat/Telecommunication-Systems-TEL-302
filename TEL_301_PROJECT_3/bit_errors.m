function num_of_bit_errors = bit_errors(est_bit_seq,b)
num_of_bit_errors=0;
for i=1:length(b)
    if(abs(b(i)-est_bit_seq(i))>0.001)
        num_of_bit_errors=num_of_bit_errors+1;
    end
end

end