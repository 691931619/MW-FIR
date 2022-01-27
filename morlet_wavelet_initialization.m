function [m1,m2, c1,c2]=morlet_wavelet_initialization(WL,fs,fc,a)
T=1/fs;
num=floor(fs*WL);
if rem(num,2)==0
    t=-num/2*T:T:num/2*T;
else
    t=-(num-1)/2*T:T:(num-1)/2*T;
end
m=exp((-(t).^2)/a);
c1=sum(m)/2;
c2=a*pi^2;
m1=exp(1i*2*pi*fc*t).*exp((-t.^2) / a) ;
m2=exp(1i*2*pi*(fc+1)*t).*exp((-t.^2) / a) ;
end