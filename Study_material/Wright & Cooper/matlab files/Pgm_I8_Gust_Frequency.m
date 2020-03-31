% Gust response of a rigid aircraft in the frequency domain – requires aircraft data, flight case and derivative codes
close all;
% RMS gust velocity (TAS) and characteristic scale wavelength (2500 ft – converted to m) 
sigma_g = 1;  L_g = 2500 / 3.2808;
 
% Equations of motion for rigid aircraft
M = [m  0; 0  I_y]; CC = [- Z_zdot  -Z_q;  -M_zdot  -M_q];  KK = [0  -Z_alpha;  0  -M_alpha];
FW = [Z_gW;  M_gW]; FT = [Z_gT;  M_gT];
 
% Calculate Gust PSD and FRF for each frequency value from 0 up to Nyquist and set Nyquist frequency 
% value to be real - note that FRF needs to be 2 dimensional (responses and frequency stored)
% The zero frequency value is NaN (Not a Number) because the KK Matrix is singular
f_nyq = 5; nt = 1024; nf = nt / 2 + 1;  frq = linspace(0, f_nyq, nf);  df = frq(2) - frq(1); 
H = zeros(2, nf); Hddot = zeros(2, nf);  i = sqrt(-1);

for ifq = 1:nf
    w = 2 * pi * frq(ifq);  omLg = w * L_g / V;  OM(ifq) = w / V;
    num = 1 + 8 / 3 * (1.339 * omLg)^2;  den = (1 + (1.339 * omLg)^2)^(11 / 6);
    psi_gr(ifq) = sigma_g^2 * L_g / pi * num / den;   		% Gust PSD for reduced frequency
    psi_gf(ifq) = sigma_g^2 * 2 * L_g / V * num / den;   		% Gust PSD for true frequency 
    HI = (KK - w^2 * M + i * w * CC)^-1;  HQG = HI * (FW + FT * exp(- w * l_WT / V));
    Hqg(:,ifq) = HQG;  Hqddotg(:,ifq) = - w^2 * Hqg(:,ifq);
end
 
% Plot gust PSD against frequency (not reduced) 
figure(1);  loglog(frq(1:nf), psi_gf(1:nf), 'k-');
xlabel('Frequency (Hz)'); ylabel('PSD of Gust Velocity    (m/s)^2/Hz');  title('Gust PSD - Frequency')

% Convert to FRF and |FRF|^2 for centre of mass response from generalised (heave / pitch) responses
Hz_C = [1  0] * Hqg;  Hz_C2 = (abs(Hz_C)).^2; Hzddot_C = [1  0] * Hqddotg;  Hzddot_C2 = (abs(Hzddot_C)).^2;
% Convert to FRF for front fuselage and tailplane response from generalised (heave / pitch) responses
Hz_F = [1  -l_F] * Hqg;  Hz_T = [1  l_T] * Hqg;
% Calculation of centre of mass response PSD from gust PSD and response-to-gust |FRF|^2  
Pzddot_C2 = psi_gf.*Hzddot_C2;
 
% Calculate root-mean-square values of acceleration at centre of mass (convert from m/s^2 to g)
g2 = 9.81^2;  sumC = 0;
for ifq = 2:nf
      sumC = sumC + Pzddot_C2(ifq) / g2 * df;
end
rmsCg = sqrt(sumC)
 
% Plot centre of mass response PSD against frequency 
figure(2);  loglog(frq(1:nf), Hzddot_C2(1:nf)/g2,'k:',frq(1:nf), Pzddot_C2(1:nf)/g2,'k-');
xlabel('Frequency (Hz)'); ylabel('FRF^2 and Acceleration PSD (g^2 / Hz)');
title('Frequency Domain Gust Response - Rigid Aircraft with Heave / Pitch Model')
legend('(Acceleration per Gust Velocity FRF)^2', 'Acceleration PSD')
