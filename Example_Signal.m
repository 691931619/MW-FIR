clc;clear;

% wavelet (filter) parameters
P_MW_FIR_WL=0.08; % P-class window length
M_MW_FIR_WL=0.28; % M-class window length
fs=800;           % sampling frequency
P_fc=42.5;        % P-class center frequency
P_fci1=112;       % enhanced P-class center frequency for interfering signal
M_fc=49.5;        % M-class center frequency 
M_fci1=25;        % enhanced M-class center frequency for interfering signal
M_fci2=35;        % enhanced M-class center frequency for interfering signal
M_fci3=65;        % enhanced M-class center frequency for interfering signal
M_fci4=75;        % enhanced M-class center frequency for interfering signal
P_a=1.55e-4;      % P-class alpha
M_a=1.2e-3;       % M-class alpha
P_lambda=0.01;    % P-class lambda
M_lambda=0.3;     % M-class lambda
R=1;              % iteration number of NSII algorithm
Q=1;              % iteration number of HIHI algorithm

% Morlet wavelet (filter) initialization
[P_m1,P_m2,P_c1,P_c2]=morlet_wavelet_initialization(P_MW_FIR_WL,fs,P_fc,P_a);
[P_mc1,P_mc2,~]=morlet_wavelet_initialization(P_MW_FIR_WL,fs,P_fci1,P_a);
[M_m1,M_m2,M_c1,M_c2]=morlet_wavelet_initialization(M_MW_FIR_WL,fs,M_fc,M_a);
[M_mc11,M_mc12,~]=morlet_wavelet_initialization(M_MW_FIR_WL,fs,M_fci1,M_a);
[M_mc21,~]=morlet_wavelet_initialization(M_MW_FIR_WL,fs,M_fci2,M_a);
[M_mc31,~]=morlet_wavelet_initialization(M_MW_FIR_WL,fs,M_fci3,M_a);
[M_mc41,M_mc42,~]=morlet_wavelet_initialization(M_MW_FIR_WL,fs,M_fci4,M_a);

% set up the input signal
amp = 1;    % fundamental amplitude
NR = 80;    % SNR, unit is dB
amp_har=0;  % interference signal amplitude 
fre_har=0;  % interference signal frequency
kx = 0;     % modulation amplitude factor
ka = 0;     % phase modulation factor
fm = 0;     % modulation frequency
pha_mod = 0;% modulation phase angle
f_offset =0;% frequency offset
Rf = 0;     % rate of change of frequency
pha = rand;    % initial phase angle
pha_har=0;  % phase angle of interference signal
endTime=50;  % signal length (needs to be changed for frequency ramping test)
t=0:1/fs:endTime;
s=amp*(1+kx*cos(2*pi*fm*t+pha_mod)).* cos(2*pi*50*t+2*pi*f_offset*t+ka*cos(2*pi*fm*t-pi+pha_mod)+pi*Rf*t.*t+pha)+amp_har*cos(2*pi*fre_har*t+pha_har);
s= awgn(s,NR,'measured');
% reference values
Tamp = amp*(1+kx*cos(2*pi*fm*t+pha_mod));
ang = 2*pi*(50+f_offset)*t + ka*cos(2*pi*fm*t-pi+pha_mod) + pi*Rf*t.*t + pha;
vector=Tamp.*exp(1i*ang);
P_TV_r=real(vector(ceil(length(P_m1)/2):1:end-floor(length(P_m1)/2)));
P_TV_i=imag(vector(ceil(length(P_m1)/2):1:end-floor(length(P_m1)/2)));
M_TV_r=real(vector(ceil(length(M_m1)/2):1:end-floor(length(M_m1)/2)));
M_TV_i=imag(vector(ceil(length(M_m1)/2):1:end-floor(length(M_m1)/2)));
fre = 50 + f_offset - fm*ka*sin(2*pi*fm*t-pi+pha_mod) + Rf*t;
P_fre=fre(ceil(length(P_m1)/2):1:end-floor(length(P_m1)/2));
M_fre=fre(ceil(length(M_m1)/2):1:end-floor(length(M_m1)/2));
ROCOF = -2*pi*fm*fm*ka*cos(2*pi*fm*t-pi+pha_mod) + Rf;
P_ROCOF=ROCOF(ceil(length(P_m1)/2)+1:1:end-floor(length(P_m1)/2)-1);
M_ROCOF=ROCOF(ceil(length(M_m1)/2)+1:1:end-floor(length(M_m1)/2)-1);

% estimation process
[P_Vector,P_Freq,P_Rocof]=MW_FIR_estimation(s,P_m1,P_m2,P_fc,P_c1,P_c2,fs);
[M_Vector,M_Freq,M_Rocof]=MW_FIR_estimation(s,M_m1,M_m2,M_fc,M_c1,M_c2,fs);
[eP_Vector,eP_Freq,eP_Rocof]=P_class_enhanced_MW_FIR_estimation(s,P_m1,P_m2,P_mc1,P_mc2,P_fc,P_fci1,P_c1,P_c2,fs,P_a,P_lambda,R,Q);
[eM_Vector,eM_Freq,eM_Rocof]=M_class_enhanced_MW_FIR_estimation(s,M_m1,M_m2,M_mc11,M_mc12,M_mc21,M_mc31,M_mc41,M_mc42,M_fc,M_fci1,M_fci2,M_fci3,M_fci4,M_c1,M_c2,fs,M_a,M_lambda,Q);

% error estimation
P_EV_r=real(P_Vector);
P_EV_i=imag(P_Vector);
P_TVE=sqrt((P_EV_r-P_TV_r).^2+(P_EV_i-P_TV_i).^2./(P_TV_r.^2+P_TV_i.^2))*100;
M_EV_r=real(M_Vector);
M_EV_i=imag(M_Vector);
M_TVE=sqrt((M_EV_r-M_TV_r).^2+(M_EV_i-M_TV_i).^2./(M_TV_r.^2+M_TV_i.^2))*100;
eP_EV_r=real(eP_Vector);
eP_EV_i=imag(eP_Vector);
eP_TVE=sqrt((eP_EV_r-P_TV_r).^2+(eP_EV_i-P_TV_i).^2./(P_TV_r.^2+P_TV_i.^2))*100;
eM_EV_r=real(eM_Vector);
eM_EV_i=imag(eM_Vector);
eM_TVE=sqrt((eM_EV_r-M_TV_r).^2+(eM_EV_i-M_TV_i).^2./(M_TV_r.^2+M_TV_i.^2))*100;
P_FE=abs(P_Freq-P_fre);
P_RFE=abs(P_Rocof-P_ROCOF);
eP_FE=abs(eP_Freq-P_fre);
eP_RFE=abs(eP_Rocof-P_ROCOF);
M_FE=abs(M_Freq-M_fre);
M_RFE=abs(M_Rocof-M_ROCOF);
eM_FE=abs(eM_Freq-M_fre);
eM_RFE=abs(eM_Rocof-M_ROCOF);

fprintf('Max.TVE for P-MW-FIR   = %d\n', max(P_TVE))
fprintf('Max.TVE for e-P-MW-FIR = %d\n', max(eP_TVE))
fprintf('Max.TVE for M-MW-FIR   = %d\n', max(M_TVE))
fprintf('Max.TVE for e-M-MW-FIR = %d\n', max(eM_TVE))

fprintf('Max.FE for P-MW-FIR    = %d\n', max(P_FE))
fprintf('Max.FE for e-P-MW-FIR  = %d\n', max(eP_FE))
fprintf('Max.FE for M-MW-FIR    = %d\n', max(M_FE))
fprintf('Max.FE for e-M-MW-FIR  = %d\n', max(eM_FE))

fprintf('Max.RFE for P-MW-FIR   = %d\n', max(P_RFE))
fprintf('Max.RFE for e-P-MW-FIR = %d\n', max(eP_RFE))
fprintf('Max.RFE for M-MW-FIR   = %d\n', max(M_RFE))
fprintf('Max.RFE for e-M-MW-FIR = %d\n', max(eM_RFE))

