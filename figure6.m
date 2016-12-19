function figure6

close all;
clear all;
clc;
%length of the signal
N=1024;
M=32;
K=128;
t=-N/2:1:N/2-1/N;
t=t/N;
T=[]; 
s=1*exp(-j*t.^2*8*M+j*t*M*16*pi)+1*exp(j*t.^2*32*M-j*t*M*16*pi)+1*exp(j*t.^2*8*M-j*t*M*16*pi);     %% 2 chirpa

x1=s.*exp(j*t.^2*18*M);     
x2=s.*exp(-j*t.^2*12*M);     
x3=s.*exp(-j*t.^2*16*M);     



for  a=[8,-32,-8] %% Change demodulation parameters in the range a=-64:1:64
xx=s.*exp(j*a*M*t.^2);   

B=dftmtx(N);
Binv=inv(B);
xf=B*xx';
q=randperm(N);
qM=q(1:K);
%creating measurement matrix
Ax=Binv(q(1:K),:);

y=(Ax*xf);
y=y';

[xp,X,Thresh]=sira(y,qM,1024,128);
T=[T Thresh];
sig_rec=ifft((xp)).*exp(-j*2*pi*a*M*t.^2);

end

xx1=zeros(size(x1));
xx2=zeros(size(x2));
xx3=zeros(size(x3));
xx1(qM)=x1(qM);
xx2(qM)=x2(qM);
xx3(qM)=x3(qM);
figure, 
%subplot(311),plot(abs(fft(s))),axis tight; 
subplot(311),plot(abs(fft(xx1)),'k'), xlabel({'Frequency samples', '(a)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(1)*ones(1,1024),'r'),axis tight;
subplot(312),plot(abs(fft(xx2)),'k'), xlabel({'Frequency samples', '(b)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(2)*ones(1,1024),'r'),axis tight;
subplot(313),plot(abs(fft(xx3)),'k'), xlabel({'Frequency samples', '(c)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(3)*ones(1,1024),'r'),axis tight;

end

function [XX1,X,Thresh]=sira(y,qM,N,M);

xr=zeros(1,N);

Amp=sum(abs(y).^2)/M;
sigma_ms2=M*(N-M)/(N-1)*Amp;
Q=1-0.99^(1/(N));

Thresh=1/M*sqrt(-(sigma_ms2)*log(Q));
X=zeros(1,N);

for k=0:N-1
    X(k+1)=1/M*sum(y.*exp(-j*2*pi*k*qM/N));
end

A=Amp;    

[a,b]=find(abs(X)>Thresh);

AFULL=dftmtx(N);
Acs=AFULL(qM,b);
X1=pinv(Acs'*Acs)*(Acs'*y');

XX1=zeros(1,N);
XX1(b)=N*X1;

end
