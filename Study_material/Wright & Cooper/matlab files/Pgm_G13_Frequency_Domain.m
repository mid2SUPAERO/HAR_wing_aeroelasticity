% Response of SDoF to single cycle of a square wave using the Fourier Transform and
% transformation into the frequency domain and back again
clear all; close all

% System parameters
f_nat = 2; period = 1 / f_nat; w_nat = 2 * pi * f_nat; zeta = 0.05;
mass = 1000; stiff = w_nat^2 * mass; force = 1000;

% Data parameters in time and frequency domains - note that time vector 
% stops dt short of T for FT analysis to remain a power of 2 and periodic
T = 8; nt = 128 ; dt = T / nt; t = [0:dt:T - dt];
df = 1 / T; f_nyq = 1 / 2 / dt; frq = [0:df:f_nyq]; [dummy,nf] = size(frq);

% Define square wave cycle excitation f(t)
pulse_width = period / 2; npulse = round(pulse_width / dt); 
for it = 1 : npulse
   f(it) = force;
end
for it = npulse + 1 : 2 * npulse
    f(it) = - force;
end
for it = 2 * npulse + 1 : nt
    f(it) = 0;
end

% Fourier Transform f(t) to frequency domain FF - apply scaling factor nt
FF = fft(f,nt) / nt;
figure(1)
plot(frq(1:nf),real(FF(1:nf)),'kx',frq(1:nf),imag(FF(1:nf)),'ko');
xlabel('Frequency (Hz)'); ylabel('Real / Imaginary Parts - FT of x(t)');
legend('Real Part','Imaginary Part')

% Generate the Frequency Response Function (FRF) for SDoF over 0 - f_nyq
for ifq = 1 : nf
    w = 2 * pi * frq(ifq);
    r = w / w_nat;
    H(ifq) = (1 / stiff)/(1 - r^2 + i * 2 * zeta * r);
end
H(nf)=real(H(nf));

% Generate the FRF 'negative' frequency content (i.e. pack to nt complex
% numbers) to be in correct format for inverse FT
for ifq = nf + 1 : nt
    H(ifq) = conj(H(nt - ifq + 2));
end
figure(2)
plot(frq(1:nf),real(H(1:nf)),'kx',frq(1:nf),imag(H(1:nf)),'ko');
xlabel('Frequency (Hz)'); ylabel('Real / Imaginary Parts - FRF');
legend('Real Part','Imaginary Part')

% Multiply FRF by FT of f(t) (element by element) to get XF - FT of x(t)
XF = H.* FF;
figure(3)
plot(frq(1:nf),real(XF(1:nf)),'kx',frq(1:nf),imag(XF(1:nf)),'ko')
xlabel('Frequency (Hz)'); ylabel('Real /Imaginary Parts - FFT of y(t)');
legend('Real Part','Imaginary Part');

% Generate response x(t) using the IFT of XF - apply scaling factor nt
x = ifft(XF) * nt;

% Plot response in mm
figure(4); subplot(211); 
plot(t,f,'k');
xlabel('Time (s)'); ylabel('Excitation (N)');
subplot(212); 
plot(t,x*1000,'k'); axis([0 T -25 25]);
xlabel('Time (s)'); ylabel('Response to Single Cycle of a Square Wave (mm)');
