clear all;
close all;
%% Question 1
N=100;
T=0.001;
A=4;
over=10;
a=0.8;
over=10;
Ts=T/over;
Nf=2048;
Fs=1/Ts;
Faxis=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;

b=(sign(randn(3*N, 1)) + 1)/2;
% T=0.001;

%% Question 2
X=bits_to_PSK_8(b);
scatterplot(X(:,1)+1i*X(:,2));
title('scatterPlot of the PSK')
length(X)
%% Question 3

[f,t]=srrc_pulse(T,over,A,a);
t2=t;
f_all_temp=zeros(length(X),length(t)+(length(X)-1)*over);
 %creating the time moved signals
 for i=0:length(X)-1
     for j=1:length(t)
         f_all_temp(i+1,j+i*over)=X(i+1,1).*f(j);
     end
 end
  %calculating the time of the signal
t1=t(end)+Ts:Ts:t(end)+T*(length(X)-1);
 t=[t t1];
 X_I=sum(f_all_temp, 1);
figure()
subplot(2,1,1)
plot(t,X_I)
title('plot of XI(t)')

hold on;


f_all_temp=zeros(length(X),length(t2)+(length(X)-1)*over);
 %creating the time moved signals
 for i=0:length(X)-1
     for j=1:length(t2)
         f_all_temp(i+1,j+i*over)=X(i+1,2)*f(j);
     end
 end
  %calculating the time of the signal
X_Q=sum(f_all_temp, 1);
subplot(2,1,2)
plot(t,X_Q)
title('plot of XQ(t)')


figure()
Px=fftshift(abs(fft(X_I,Nf)).^2)*Ts./length(t);
subplot(2,1,1)
semilogy(Faxis,Px)
title('Plot of the periodogram of X_I(t)')
subplot(2,1,2)
Px=fftshift(abs(fft(X_Q,Nf)).^2)*Ts./length(t);
semilogy(Faxis,Px)
title('Plot of the periodogram of X_Q(t)')


%% Question 4
F0=2000;

X_I_new=(X_I.*2.*cos(2.*pi.*F0.*t));
X_Q_new=(X_Q.*(-2*sin(2.*pi.*F0.*t)));
figure()
subplot(2,1,1)
plot(t,X_I_new)
title('plot of XI(t)')

subplot(2,1,2)
plot(t,X_Q_new)
title('plot of XQ(t)')

figure()
Px=fftshift(abs(fft(X_I_new,Nf)).^2)*Ts./length(t);
subplot(2,1,1)
semilogy(Faxis,Px)
title('Periodogram of XI(t)')

subplot(2,1,2)
Px=fftshift(abs(fft(X_Q_new,Nf)).^2)*Ts./length(t);
semilogy(Faxis,Px)
title('Periodogram of XQ(t)')

%% Question 5

Xt=X_I_new+X_Q_new;
figure()
subplot(2,1,1)
plot(t, Xt)
title('plot of X(t)')

Px=fftshift(abs(fft(Xt,Nf)).^2)*Ts./length(t);
subplot(2,1,2)
semilogy(Faxis,Px)
title('Periodogram of XI(t)')

%% Question 7

SNRdb=10;
varw=(1/Ts.*10^(SNRdb/10));
varN=Ts*varw/2;
Gaussian_Noise=sqrt(varN).*randn(1,length(Xt));

Yt=Xt+Gaussian_Noise;
figure()
plot(t,Yt)
title('plot of X(t) with Noise')

%% Question 8

Y_I_rec = Yt.*cos(2*pi*F0*t);
Y_Q_rec = -Yt.*sin(2*pi*F0*t);

figure()
subplot(2,1,1)
plot(t, Y_I_rec)
title('plot of YI(t)')

subplot(2,1,2)
plot(t,Y_Q_rec)
title('plot of YQ(t)')


Px=fftshift(abs(fft(Y_I_rec,Nf)).^2)*Ts./length(t);
figure()
subplot(2,1,1)
semilogy(Faxis, Px)
title('Periodogram of YI(t)')

Px=fftshift(abs(fft(Y_Q_rec,Nf)).^2)*Ts./length(t);
subplot(2,1,2)
semilogy(Faxis, Px)
title('Periodogram of YQ(t)')

%% Question 9

Y_I_filtered=conv(f,Y_I_rec)*Ts;
Y_Q_filtered=conv(f,Y_Q_rec)*Ts;

t_conv=t(1)+t2(1):Ts:t(end)+t2(end);
figure()
subplot(2,1,1)
plot(t_conv,Y_I_filtered)
title('Plot of Y_I(t) filtered')

subplot(2,1,2)
plot(t_conv,Y_Q_filtered)
title('Plot of Y_Q(t) filtered')

Px=fftshift(abs(fft(Y_I_filtered,Nf)).^2)*Ts./length(t);
figure()
subplot(2,1,1)
semilogy(Faxis, Px)
title('Periodogram of YI(t) filtered')

Px=fftshift(abs(fft(Y_Q_filtered,Nf)).^2)*Ts./length(t);
subplot(2,1,2)
semilogy(Faxis, Px)
title('Periodogram of YQ(t) filtered')

%% Question 10


j=0;
while t_conv(j+1)<0
     j=j+1;  
end
 t_temp=zeros(1,N);
 Y_I_final=zeros(1,N);
 Y_Q_final=zeros(1,N);

 
for i=1:N
    t_temp(i)=t_conv(j+(i-1)*over);
    Y_I_final(i)=Y_I_filtered(j+(i-1)*over);
    Y_Q_final(i)=Y_Q_filtered(j+(i-1)*over);
end


scatterplot(Y_I_final+1i*Y_Q_final);
title('The YI and QI final output')

%% Question 11

Y_Final(:,1)=Y_I_final(1:N);
Y_Final(:,2)=Y_Q_final(1:N);
[est_X,est_bit_seq]=detect_PSK_8(Y_Final);
%% Question 12

symbol_errors(est_X,X)

%% Question 13

bit_errors(est_bit_seq,b)

%% B meros

SNR=[-2:2:16];

k=200;
Esymbol=zeros(length(SNR),1);
Ebit=zeros(length(SNR),1);

 
for p=1:10
       sum1=0;
       sum2=0;
   for  u=1:k
    b=(sign(randn(3*N, 1)) + 1)/2;
    
    X=bits_to_PSK_8(b);
    [f,t]=srrc_pulse(T,over,A,a);
    t2=t;
    f_all_temp=zeros(length(X),length(t)+(length(X)-1)*over);
 	%creating the time moved signals
     for i=0:length(X)-1
         for j=1:length(t)
         	 f_all_temp(i+1,j+i*over)=X(i+1,1).*f(j);
         end
     end
    %calculating the time of the signal
    t1=t(end)+Ts:Ts:t(end)+T*(length(X)-1);
    t=[t t1];
    X_I=sum(f_all_temp, 1);
    
    f_all_temp=zeros(length(X),length(t2)+(length(X)-1)*over);
    %creating the time moved signals
     for i=0:length(X)-1
          for j=1:length(t2)
             f_all_temp(i+1,j+i*over)=X(i+1,2)*f(j);
          end
     end
         %calculating the time of the signal
    X_Q=sum(f_all_temp, 1);
    
    
    X_I_new=(X_I.*2.*cos(2.*pi.*F0.*t));
    X_Q_new=(X_Q.*(-2*sin(2.*pi.*F0.*t)));  
    
    Xt=X_I_new+X_Q_new;
    
    SNRdb=SNR(p);
    
    
    varw=(1/Ts.*10^(SNRdb/10));
    varN=Ts*varw/2;
    Gaussian_Noise=sqrt(varN).*randn(1,length(Xt));
    Yt=Xt+Gaussian_Noise;
    
    Y_I_rec = Yt.*cos(2*pi*F0*t);
    Y_Q_rec = -Yt.*sin(2*pi*F0*t);
    
    Y_I_filtered=conv(f,Y_I_rec)*Ts;
    Y_Q_filtered=conv(f,Y_Q_rec)*Ts;

    t_conv=t(1)+t2(1):Ts:t(end)+t2(end);    
    
    j=0;
    while t_conv(j+1)<0
             j=j+1;  
    end
    t_temp=zeros(1,N);
    Y_I_final=zeros(1,N);
    Y_Q_final=zeros(1,N);

 
    for i=1:N
     t_temp(i)=t_conv(j+(i-1)*over);
     Y_I_final(i)=Y_I_filtered(j+(i-1)*over);
     Y_Q_final(i)=Y_Q_filtered(j+(i-1)*over);
    end
    Y_Final(:,1)=Y_I_final(1:N);
    Y_Final(:,2)=Y_Q_final(1:N);
    [est_X,est_bit_seq]=detect_PSK_8(Y_Final);

    sum1=sum1+symbol_errors(est_X,X);

    sum2=sum2+bit_errors(est_bit_seq,b);
  

   end
    Esymbol(p)=sum1/(N*k);
    Ebit(p)=sum2/(N*3*k);
end

temp=Q(SNR);
figure()
semilogy(SNR, Esymbol)
hold on
semilogy(SNR, temp)
