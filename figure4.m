function figure4

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
s=1*exp(-j*t.^2*8*M+j*t*M*16*pi)+1*exp(j*t.^2*32*M-j*t*M*16*pi)+1*exp(j*t.^2*8*M-j*t*M*16*pi);     %% 2 chirps

x1=s.*exp(j*t.^2*8*M);     
x2=s.*exp(-j*t.^2*32*M);     
x3=s.*exp(-j*t.^2*8*M);     


figure(1), %% Fig. 2 
subplot(411),plot(abs(fft(s)),'k'),axis tight; xlabel({'Frequency samples', '(a)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(412),plot(abs(fft(x1)),'k'),axis tight; xlabel({'Frequency samples', '(b)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(413),plot(abs(fft(x2)),'k'),axis tight; xlabel({'Frequency samples', '(c)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(414),plot(abs(fft(x3)),'k'),axis tight; xlabel({'Frequency samples', '(d)'}), ylabel('Amplitude'), colormap(1-gray) 

end


