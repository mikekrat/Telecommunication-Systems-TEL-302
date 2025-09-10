function X=bits_to_4PAM(c,d)


for j=2:2:length(c)
    if c(j-1)==0 && c(j)==0
        d(j/2)=3;
    elseif c(j-1)==0 && c(j)==1
        d(j/2)=1;
    elseif c(j-1)==1 && c(j)==1
        d(j/2)=-1;     
    else
        d(j/2)=-3;
        
    end
end

X=d;
end


