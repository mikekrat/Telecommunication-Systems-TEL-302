clear all;
close all;
%% A.1
T=1/1000;
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=[0 0.5 1];
Nf=2048;
%we create 3 SRRC pulses
[f1,t] = srrc_pulse(T,over,A,a(1));
[f2,t] = srrc_pulse(T,over,A,a(2));
[f3,t] = srrc_pulse(T,over,A,a(3));

%we plot them on the same figure
figure()
plot(t,f1,'r')
title('Plot of SRRC Signals')

hold on;
plot(t,f2,'green')
hold on;
plot(t,f3,'yellow')
%% A.2
%Calculating FFT and it's time axis
Faxis=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
F1=fftshift(fft(f1,Nf))*Ts;
F2=fftshift(fft(f2,Nf))*Ts;
F3=fftshift(fft(f3,Nf))*Ts;
figure()
plot(Faxis,abs(F1).^2)
title('Energy Spectral Density');
hold on;
plot(Faxis,abs(F2).^2)
hold on;
plot(Faxis,abs(F3).^2)
figure()
semilogy(Faxis,abs(F1).^2)
title('Energy Spectral Density');
hold on;
semilogy(Faxis,abs(F2).^2)
hold on;
semilogy(Faxis,abs(F3).^2)
hold on;

%% A.3

BW=(1+a)/(2*T);
disp('The Theoretical BW of F1 is: ');
disp(BW(1));
disp('The Theoretical BW of F2 is: ');
disp(BW(2));
disp('The Theoretical BW of F3 is: ');
disp(BW(3));

c=T/10^3;
plot(xlim, [1 1]*c, '--k')   
hold on;
c=T/10^5;
plot(xlim, [1 1]*c, '--k')  
hold off;





%initialization
T=1/1000;
over=10;
Ts=T/over;
A=5;
a=[0 0.5 1];
%pre allocating space for the integral calculation
integral=zeros(length(a),4);
k=0:2*A;

%% B.1.1
%for each value of a
for u=1:length(a)
 tnew=[-A*T:Ts:A*T + k(4)*T] + 10^(-8);    
[f,t]=srrc_pulse(T,over,A,a(u));
ftemp=zeros(1,length(tnew)-length(t));
%The vector is changed to the length of the time shifted f's
f=[f ftemp];
  figure()
    plot(tnew,f)
    hold on;
    f_all=zeros(4,length(tnew));
for i=k(1):k(4)
    %creating the time axis for each time shifted signal
    tnew=[-A*T:Ts:A*T + i*T] + 10^(-8);
    for j=1 : length(t)
        offset=j+length(tnew)-length(t);
       f_all(i+1,offset)=f(j); 
    end
    %pre allocating space for each time shifted signal
    hey=zeros(1,length(tnew));
    for j=1 : length(hey)
        hey(j)=f_all(i+1,j); 
    end
    %plotting the time shifted signals
    plot(tnew,hey);
    hold on;
end
hold off;

%% B.1.2
figure()
%claculating f(t)*f(t-kT) for each k
for i=k(1):k(4)
    hey=zeros(1,length(tnew));
    for j=1 : length(hey)
        hey(j)=f_all(i+1,j); 
    end 
    answer=hey.*f;
    plot(tnew,answer)
    hold on;
    %calculating the integrals
    answer=answer.*Ts;
    integral(u,i+1)=sum(answer);
end
hold off;
end
  %% B.1.3
disp('THe value of the integrals are:');
disp('For k=0: ')
disp(integral(:,1));
disp('For k=1: ')
disp(integral(:,2));
disp('For k=2: ')
disp(integral(:,3));
disp('For k=3: ')
disp(integral(:,4));




%% C.1

T=0.1;
over=10;
Ts=T/over;
a=0.5;
A=5;
N=50;

[f2,t]=srrc_pulse(T,over,A,a);
b = (sign(randn(N, 1)) + 1)/2;
%% C.2 
X = bits_to_2PAM(b);
X_delta =1/Ts*upsample(X,over);
n=0:Ts:N*T-Ts;
figure()
stem(n,X_delta*Ts)
figure()
%calculating the convolution and it's time axis
Xconv=conv(f2,X_delta)*Ts;
t_Xconv=linspace(t(1)+n(1),t(end)+n(end),length(Xconv));
plot(t_Xconv,Xconv)

%creating f(-t)
f2rever=fliplr(f2);
f2rever_time=-fliplr(t);
%calculating the convolution and it's time axis
Z=conv(Xconv,f2rever)*Ts;
t_Z=linspace(f2rever(1)+t_Xconv(1),f2rever(end)+t_Xconv(end),length(Z));
figure()
plot(t_Z,Z)
hold on;
temp=[0 : N - 1] * T;
stem(temp, X)

