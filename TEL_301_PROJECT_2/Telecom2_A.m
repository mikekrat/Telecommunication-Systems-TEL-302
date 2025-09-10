clear all;
close all;
%% A
%value initialization
T=1/100;
over=10;
Ts=T/over;
A=4;
a=0.5;
Nf=2048;
Fs=1/Ts;
N=100;

%% A.1
%creating an srrc pulse and calculating it's power sprectral density
[f1,t1]=srrc_pulse(T,over,A,a);
Faxis=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
F1=abs(fftshift(fft(f1,Nf))*Ts);
figure()
semilogy(Faxis,F1.^2)
title('Plot of the FFt of the srrc pulse ')
%% A.2
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
 %calculating the time of the signal
 t=t1(end)+Ts:Ts:t1(end)+T*(length(b)-1);
 t=[t1 t];
 X_t=sum(f_all_temp, 1);
 
 figure()
 plot(t,X_t)
title('Plot of X(t) constructed with 2-PAM ')
%% A.3
figure()

%creating the periodogram of X_t
Px=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
subplot(2,1,1)
plot(Faxis,Px)
title('Plot of the periodogram of X(t)(standard plot)')
subplot(2,1,2)
semilogy(Faxis,Px)

title('Plot of the periodogram of X(t)(semilogy)')
%calculating the PSD of Px
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
 Pmean(u,:)=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
 
end
%summing all the calculated X_t and calculating the practical PSD of X_t
PSD_pract=sum(Pmean,1)/rep;
%Calculating the theoretical PSD of X_t
PSD_theor=(1/T).*(F1.^2);
figure()
semilogy(Faxis,PSD_pract);
hold on;
semilogy(Faxis,PSD_theor);
title('Plot of the theoretical and practical PSD of X(t)')
%% A.4
%creating an array of 100 random bits
b = (sign(randn(N, 1)) + 1)/2;
d = zeros(1,length(b)/2);
X = bits_to_4PAM(b,d);
 %creating an empty vector to store the numbers
 f_all_temp=zeros(length(X),length(t1)+(length(X)-1)*over);
 %creating the time moved signals
 for i=0:length(X)-1
     for j=1:length(t1)
         f_all_temp(i+1,j+i*over)=X(i+1)*f1(j);
     end
 end
  %calculating the time of the signal
t=t1(end)+Ts:Ts:t1(end)+T*(length(X)-1);
 t=[t1 t];
 X_t=sum(f_all_temp, 1);
figure()
Px=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
subplot(2,1,1)
plot(Faxis,Px)
title('Plot of the periodogram of X(t)(standard plot)')
subplot(2,1,2)
semilogy(Faxis,Px)
title('Plot of the periodogram of X(t)(semilogy)')
rep=500;
Pmean=zeros(rep,length(Px));

for u=1:rep
%creating an array of 100 random bits
b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);
 %creating an empty vector to store the numbers
 f_all_temp=zeros(length(X),length(t1)+(length(X)-1)*over);
 %creating the time moved signals
 for i=0:length(X)-1
     for j=1:length(t1)
         f_all_temp(i+1,j+i*over)=X(i+1)*f1(j);
     end
 end
 
 X_t=sum(f_all_temp, 1);
 Pmean(u,:)=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
 
end
%summing all the calculated X_t and calculating the practical PSD of X_t
PSD_pract=sum(Pmean,1)/rep;
%Calculating the theoretical PSD of X_t
PSD_theor=(1/T).*(F1.^2);
figure()
semilogy(Faxis,PSD_pract);

hold on;
semilogy(Faxis,PSD_theor);
title('Plot of the theoretical and practical PSD of X(t)')
%% A.5

T=2*T;
over=20;
Ts=T/over;
[f1,t1]=srrc_pulse(T,over,A,a);
Ts=T/over;
Fs=1/Ts;
Faxis=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
F1=abs(fftshift(fft(f1,Nf))*Ts);



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

figure()
 plot(t,X_t)
 title('Plot of X(t) constructed with 2-PAM ')
 
 
figure()
%creating the periodogram of X_t
Px=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
subplot(2,1,1)
plot(Faxis,Px)
title('Plot of the periodogram of X(t)(standard plot) for T=2T')
subplot(2,1,2)
semilogy(Faxis,Px)
title('Plot of the periodogram of X(t)(semilogy) for T=2T')
 
 
 
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
 Pmean(u,:)=fftshift(abs(fft(X_t,Nf)).^2)*Ts./length(t);
end
%summing all the calculated X_t and calculating the practical PSD of X_t
PSD_pract=sum(Pmean,1)/rep;
%Calculating the theoretical PSD of X_t
PSD_theor=(1/T).*(F1.^2);

figure()
semilogy(Faxis,PSD_pract);
hold on;
semilogy(Faxis,PSD_theor);
title('Plot of the theoretical and practical PSD of X(t)')

