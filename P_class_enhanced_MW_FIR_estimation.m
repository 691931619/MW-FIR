function [Vector,Freq,Rocof]=P_class_enhanced_MW_FIR_estimation(s,m1,m2,mc1,mc2,fc,fcc,c1,c2,fs,a,lambda,R,Q)
% convolution between the input signal and the wavelets (filters)
cy1=conv(s,m1,'valid');
cy2=conv(s,m2,'valid');
cy1t=zeros(1,length(cy1));
cy2t=zeros(1,length(cy2));
cyc1=conv(s,mc1,'valid');
cyc2=conv(s,mc2,'valid');
cyc2t=zeros(1,length(cyc2));

% conventional MW-FIR estimation
ro=log(abs(cy1)./abs(cy2));
Freq=fc+1/2-ro/(2*c2);
phase1=angle(cy1);
phase2=angle(cy2);
ang1=phase1;
amplitude=abs(cy1)./(c1*exp((-c2)*(Freq-fc).^2));
Vector=amplitude.*exp(1i*phase1);

% remove the negative image interference
for iteration=1:1:R    
    neg_amp1=amplitude*c1.*exp(-a*pi^2*(Freq+fc).^2);
    neg_amp2=amplitude*c1.*exp(-a*pi^2*(Freq+fc+1).^2);
    ang1=phase1;
    ang2=phase2;
    neg_vec1=neg_amp1.*conj(exp(1i*ang1));
    neg_vec2=neg_amp2.*conj(exp(1i*ang2));
    cy1f=cy1-neg_vec1;
    cy2f=cy2-neg_vec2;
    ro=log(abs(cy1f)./abs(cy2f));
    Freq=fc+1/2-ro/(2*c2);
    phase1=angle(cy1f);
    phase2=angle(cy2f);
    amplitude=abs(cy1f)./(c1*exp((-c2)*(Freq-fc).^2));
    Vector=amplitude.*exp(1i*phase1);
end

% compute the ratio E
vector=amplitude*c1.*exp(-c2*(Freq-fcc).^2).*(exp(1i*ang1))+amplitude*c1.*exp(-c2*(Freq+fcc).^2).*conj(exp(1i*ang1));
har_c=abs(cyc1-vector).^2;
cyc1t=har_c./abs(cyc1);
index=find(cyc1t>=lambda);

% find the location where ratio E exceeds the threshold parameter lambda
% and update the estimates by eliminating the interference signal
for i=1:length(index)
    for o=1:1:Q
        ti=index(i);
        compare_bin1=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc).^2).*(exp(1i*phase1(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc).^2).*conj(exp(1i*phase1(ti)));
        compare_bin2=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc-1).^2).*(exp(1i*phase2(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc+1).^2).*conj(exp(1i*phase2(ti)));
        cyc1t(ti)=cyc1(ti)-compare_bin1;
        cyc2t(ti)=cyc2(ti)-compare_bin2;
        ro_compare=log(abs(cyc1t(ti))./abs(cyc2t(ti)));
        Freq_compare=fcc+1/2-ro_compare/(2*c2);
        phase1_compare=angle(cyc1t(ti));
        phase2_compare=angle(cyc2t(ti));    
        amplitude_compare=abs(cyc1t(ti))./(c1*exp((-c2)*(Freq_compare-fcc).^2));

        vector_compare1=amplitude_compare*c1.*exp(-a*pi^2*(Freq_compare-fc).^2).*(exp(1i*phase1_compare))+amplitude_compare*c1.*exp(-a*pi^2*(Freq_compare+fc).^2).*conj(exp(1i*phase1_compare));
        vector_compare2=amplitude_compare*c1.*exp(-a*pi^2*(Freq_compare-fc-1).^2).*(exp(1i*phase2_compare))+amplitude_compare*c1.*exp(-a*pi^2*(Freq_compare+fc+1).^2).*conj(exp(1i*phase2_compare));
        cy1t(ti)=cy1(ti)-vector_compare1;
        cy2t(ti)=cy2(ti)-vector_compare2;

        ro=log(abs(cy1t(ti))./abs(cy2t(ti)));
        Freq(ti)=fc+1/2-ro/(2*c2);
        phase1(ti)=angle(cy1t(ti));
        amplitude(ti)=abs(cy1t(ti))./(c1*exp((-c2)*(Freq(ti)-fc).^2));
        Vector(ti)=amplitude(ti).*exp(1i*phase1(ti));
    end
end
% compute the ROCOF
Rocof_1=Freq(1:end-2);
Rocof_2=Freq(3:end);
Rocof=(Rocof_2-Rocof_1)/2/(1/fs);
end