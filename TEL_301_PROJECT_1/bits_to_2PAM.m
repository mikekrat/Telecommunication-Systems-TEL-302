function X=bits_to_2PAM(b)

for i=1:length(b)
    if b(i)==0
        b(i)=1;
    else
        b(i)=-1;
    end
end

X=b;
end

