function [Vector,Freq,Rocof]=MW_FIR_estimation(s,m1,m2,fc,c1,c2,fs)
cy1=conv(s,m1,'valid');
cy2=conv(s,m2,'valid');

ro=log(abs(cy1)./abs(cy2));
Freq=fc+1/2-ro/(2*c2);
phase=angle(cy1);
amplitude=abs(cy1)./(c1*exp((-c2)*(Freq-fc).^2));
Vector=amplitude.*exp(1i*phase);
Rocof_1=Freq(1:end-2);
Rocof_2=Freq(3:end);
Rocof=(Rocof_2-Rocof_1)/2/(1/fs);

end