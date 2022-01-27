function [Vector,Freq,Rocof]=M_class_enhanced_MW_FIR_estimation(s,m1,m2,mc11,mc12,mc21,mc31,mc41,mc42,fc,fcc1,fcc2,fcc3,fcc4,c1,c2,fs,a,lambda,Q)
% convolution between the input signal and the wavelets (filters)
cy1=conv(s,m1,'valid');
cy2=conv(s,m2,'valid');
cy1t=zeros(1,length(cy1));
cy2t=zeros(1,length(cy2));
cyc11=conv(s,mc11,'valid');
cyc12=conv(s,mc12,'valid');
cyc11t=zeros(1,length(cyc11));
cyc12t=zeros(1,length(cyc12));
cyc21=conv(s,mc21,'valid');
cyc31=conv(s,mc31,'valid');
cyc41=conv(s,mc41,'valid');
cyc42=conv(s,mc42,'valid');
cyc41t=zeros(1,length(cyc41));
cyc42t=zeros(1,length(cyc42));

% conventional MW-FIR estimation
ro=log(abs(cy1)./abs(cy2));
Freq=fc+1/2-ro/(2*c2);
phase1=angle(cy1);
ang1=phase1;
phase2=angle(cy2);
amplitude=abs(cy1)./(c1*exp((-c2)*(Freq-fc).^2));
Vector=amplitude.*exp(1i*phase1);

% compute the ratio E
vector=amplitude*c1.*exp(-c2*(Freq-fcc1).^2).*(exp(1i*ang1))+amplitude*c1.*exp(-c2*(Freq+fcc1).^2).*conj(exp(1i*ang1));
har_c=abs(cyc11-vector).^2;
cyc11_temp=har_c./abs(cyc11);
vector=amplitude*c1.*exp(-c2*(Freq-fcc2).^2).*(exp(1i*ang1))+amplitude*c1.*exp(-c2*(Freq+fcc2).^2).*conj(exp(1i*ang1));
har_c=abs(cyc21-vector).^2;
cyc21_temp=har_c./abs(cyc21);
vector=amplitude*c1.*exp(-c2*(Freq-fcc3).^2).*(exp(1i*ang1))+amplitude*c1.*exp(-c2*(Freq+fcc3).^2).*conj(exp(1i*ang1));
har_c=abs(cyc31-vector).^2;
cyc31_temp=har_c./abs(cyc31);
vector=amplitude*c1.*exp(-c2*(Freq-fcc4).^2).*(exp(1i*ang1))+amplitude*c1.*exp(-c2*(Freq+fcc4).^2).*conj(exp(1i*ang1));
har_c=abs(cyc41-vector).^2;
cyc41_temp=har_c./abs(cyc41);

% find the location where ratio E exceeds the threshold parameter lambda
% and update the estimates by eliminating the interference signal
index1=find(cyc11_temp>=lambda);
index2=find(cyc21_temp>=lambda);
index3=find(cyc31_temp>=lambda);
index4=find(cyc41_temp>=lambda);
[~,n]=max([cyc11_temp;cyc21_temp;cyc31_temp;cyc41_temp],[],1);
for i=1:length(cy1)
    ti=i;    
    if (ismember(i,index1)==1 || ismember(i,index2)==1) && (n(i)==1 || n(i)==2)
        for o=1:1:Q
            compare_bin1=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc1).^2).*(exp(1i*phase1(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc1).^2).*conj(exp(1i*phase1(ti)));
            compare_bin2=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc1-1).^2).*(exp(1i*phase2(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc1+1).^2).*conj(exp(1i*phase2(ti)));
            cyc11t(ti)=cyc11(ti)-compare_bin1;
            cyc12t(ti)=cyc12(ti)-compare_bin2;
            ro_compare=log(abs(cyc11t(ti))./abs(cyc12t(ti)));
            Freq_compare=fcc1+1/2-ro_compare/(2*c2);
            phase1_compare=angle(cyc11t(ti));
            phase2_compare=angle(cyc12t(ti));    
            amplitude_compare=abs(cyc11t(ti))./(c1*exp((-c2)*(Freq_compare-fcc1).^2));

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
    elseif (ismember(i,index3)==1 || ismember(i,index4)==1) && (n(i)==3 || n(i)==4)
        for o=1:1:Q
            compare_bin1=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc4).^2).*(exp(1i*phase1(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc4).^2).*conj(exp(1i*phase1(ti)));
            compare_bin2=amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)-fcc4-1).^2).*(exp(1i*phase2(ti)))+amplitude(ti)*c1.*exp(-a*pi^2*(Freq(ti)+fcc4+1).^2).*conj(exp(1i*phase2(ti)));
            cyc41t(ti)=cyc41(ti)-compare_bin1;
            cyc42t(ti)=cyc42(ti)-compare_bin2;
            ro_compare=log(abs(cyc41t(ti))./abs(cyc42t(ti)));
            Freq_compare=fcc4+1/2-ro_compare/(2*c2);
            phase1_compare=angle(cyc41t(ti));
            phase2_compare=angle(cyc42t(ti));    
            amplitude_compare=abs(cyc41t(ti))./(c1*exp((-c2)*(Freq_compare-fcc4).^2));

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
end
% compute the ROCOF
Rocof_1=Freq(1:end-2);
Rocof_2=Freq(3:end);
Rocof=(Rocof_2-Rocof_1)/2/(1/fs);
end