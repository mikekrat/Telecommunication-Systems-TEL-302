clear all;
close all;
T=1/100;
over=10;
Ts=T/over;
A=4;
a=0.5;
Nf=2048;
Fs=1/Ts;
N=100;
fo=1/T;


[f1,t1]=srrc_pulse(T,over,A,a);
Faxis=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
F1=abs(fftshift(fft(f1,Nf))*Ts);
theta=0+(2*pi)*rand(1,1);


b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);
%creating an empty vector to store the numbers
 f_all_temp=zeros(length(b),length(t1)+(length(b)-1)*over);
 %creating the time moved signals
 for i=0:length(b)-1
     for j=1:length(t1)
         f_all_temp(i+1,j+i*over)=X(i+1)*f1(j);
     end
 end
 %calculating the time of the signal
 t=t1(end)+Ts:Ts:t1(end)+T*(length(b)-1);
 t=[t1 t];
 X_t=sum(f_all_temp, 1);
 temp=cos(2*pi*fo*t+theta);
 Y_t=X_t.*temp;
 figure()
 plot(t,Y_t)
title('Plot of Y(t) constructed with 2-PAM ')
%creating the periodogram of Y_t
Px=fftshift(abs(fft(Y_t,Nf)).^2)*Ts./length(t);
rep=500;
Pmean=zeros(rep,length(Px));

for u=1:rep
%creating an array of 100 random bits
b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);
%creating an empty vector to store the numbers
 f_all_temp=zeros(length(b),length(t1)+(length(b)-1)*over);
 %creating the time moved signals
 for i=0:length(b)-1
     for j=1:length(t1)
         f_all_temp(i+1,j+i*over)=X(i+1)*f1(j);
     end
 end
 
 X_t=sum(f_all_temp, 1);
 temp=cos(2*pi*fo*t+theta);
 Y_t=X_t.*temp;
 Pmean(u,:)=fftshift(abs(fft(Y_t,Nf)).^2)*Ts./length(t);
 
end
%summing all the calculated Y_t and calculating the practical PSD of Y_t
PSD_pract=sum(Pmean,1)/rep;
%Calculating the theoretical PSD of Y_t
PSD_theor=(1/T).*(F1.^2);
figure()
semilogy(Faxis,PSD_pract);

title('Plot of the theoretical and practical PSD of Y(t)')
