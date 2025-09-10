function [est_X,est_bit_seq]=detect_PSK_8(Y)

est_X=zeros(length(Y),1);
est_bit_seq=zeros(3*length(Y),1);
for k=1:length(Y)

    if(Y(k,1)>(cosd(337.5)))%0 degrees
        est_X(k,1)=1;
        est_X(k,2)=0;
        est_bit_seq(k*3-2)=0;
        est_bit_seq(k*3-1)=0;
        est_bit_seq(k*3)=0;
        
    elseif(cosd(67.5)<Y(k,1)&&Y(k,1)<cosd(22.5)&&Y(k,2)>0)%45 degrees
        est_X(k,1)=1/sqrt(2);
        est_X(k,2)=1/sqrt(2);
        est_bit_seq(k*3-2)=0;
        est_bit_seq(k*3-1)=0;
        est_bit_seq(k*3)=1;
        
    elseif(Y(k,2)>sind(67.5))%90 degrees
        est_X(k,1)=0;
        est_X(k,2)=1;
        est_bit_seq(k*3-2)=0;
        est_bit_seq(k*3-1)=1;
        est_bit_seq(k*3)=1;

        
        
    elseif(cosd(112.5)>Y(k,1)&& Y(k,1)>cosd(157.5)&&Y(k,2)>0)%135 degrees
        est_X(k,1)=-1/sqrt(2);
        est_X(k,2)=1/sqrt(2);
        est_bit_seq(k*3-2)=0;
        est_bit_seq(k*3-1)=1;
        est_bit_seq(k*3)=0;

        
        
        
    elseif(Y(k,1)<cosd(157.5))%180 degrees 
        est_X(k,1)=-1;
        est_X(k,2)=0;
        est_bit_seq(k*3-2)=1;
        est_bit_seq(k*3-1)=1;
        est_bit_seq(k*3)=0;
        
        
    elseif(cosd(202.5)<Y(k,1)&& Y(k,1)<cosd(247.5) && Y(k,2)<0)%225 degrees
        est_X(k,1)=-1/sqrt(2);
        est_X(k,2)=-1/sqrt(2);
        est_bit_seq(k*3-2)=1;
        est_bit_seq(k*3-1)=1;
        est_bit_seq(k*3)=1;
        
        
    elseif(Y(k,2)<sind(247.5))%270 degrees
        est_X(k,1)=0;
        est_X(k,2)=-1;
        est_bit_seq(k*3-2)=1;
        est_bit_seq(k*3-1)=0;
        est_bit_seq(k*3)=1;
        
        
        
        %est_X(k,1)>cos(292.5)&&est_X(k,1)<cos(337,5)&&est_X(k,2)<0
    elseif(Y(k,1)>cosd(292.5) && Y(k,1)<cosd(337.5) && Y(k,2)<0)
        est_X(k,1)=1/sqrt(2);
        est_X(k,2)=-1/sqrt(2);
        est_bit_seq(k*3-2)=1;
        est_bit_seq(k*3-1)=0;
        est_bit_seq(k*3)=0;
    end
end