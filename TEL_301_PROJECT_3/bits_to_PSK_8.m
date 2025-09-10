function X=bits_to_PSK_8(bit_seq)

N=length(bit_seq)/3;
X=zeros(N,2);

for k = 0:3:length(bit_seq)-1 
    i=k/3+1;
    if((bit_seq(k+1)==0) && (bit_seq(k+2)==0) && bit_seq(k+3)==0)
    X(i,1)=cos(0);
    X(i,2)=sin(0);
    elseif((bit_seq(k+1)==0) && (bit_seq(k+2)==0) && bit_seq(k+3)==1)
    X(i,1)=cos(2*pi*1/8);
    X(i,2)=sin(2*pi*1/8);   
     elseif((bit_seq(k+1)==0) && (bit_seq(k+2)==1) && bit_seq(k+3)==1)
    X(i,1)=cos(2*pi*2/8);
    X(i,2)=sin(2*pi*2/8);    
     elseif((bit_seq(k+1)==0) && (bit_seq(k+2)==1) && bit_seq(k+3)==0)
    X(i,1)=cos(2*pi*3/8);
    X(i,2)=sin(2*pi*3/8);   
     elseif((bit_seq(k+1)==1) && (bit_seq(k+2)==1) && bit_seq(k+3)==0)
    X(i,1)=cos(2*pi*4/8);
    X(i,2)=sin(2*pi*4/8);   
     elseif((bit_seq(k+1)==1) && (bit_seq(k+2)==1) && bit_seq(k+3)==1)
    X(i,1)=cos(2*pi*5/8);
    X(i,2)=sin(2*pi*5/8);    
     elseif((bit_seq(k+1)==1) && (bit_seq(k+2)==0) && bit_seq(k+3)==1)
    X(i,1)=cos(2*pi*6/8);
    X(i,2)=sin(2*pi*6/8) ;   
     elseif((bit_seq(k+1)==1) && (bit_seq(k+2)==0) && bit_seq(k+3)==0)
    X(i,1)=cos(2*pi*7/8);
    X(i,2)=sin(2*pi*7/8) ;   
    end
end

end