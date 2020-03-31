% Response of SDoF system to single cycle of a square wave
clear all; close all;

% System parameters
f_nat = 2; period = 1 / f_nat; w_nat = 2 * pi * f_nat; zeta = 0.05;
mass = 1000; stiff = w_nat^2 * mass; force = 1000;
w_dpd = w_nat * sqrt(1 - zeta^2); psi = atan2(sqrt(1 - zeta^2) , zeta);

% Data parameters
T = 8; dt = 0.01; t = [0:dt:T]; [dummy,nt] = size(t);

% Response to step force 
for it = 1:nt
    a = exp(-zeta * w_nat * t(it)) / sqrt(1 - zeta^2);
    b = sin(w_dpd * t(it) + psi);
    s(it)= force / stiff * (1 - a * b);
end

% Response to square wave using superposition
% Function p shows 'shape' of excitation force
pulse_width = period / 2; npulse = round(pulse_width / dt); 
for it = 1 : npulse
    x(it) = s(it); f(it) = 10;
end
for it = npulse + 1 : 2 * npulse
    x(it) = s(it) - 2 * s(it - npulse);
    f(it) = - 10;
end
for it = 2 * npulse + 1 : nt
    x(it) = s(it) - 2 * s(it - npulse) + s(it - 2 * npulse);
    f(it) = 0;
end

% Plot response in mm 
plot(t,x*1000,'k-',t,f,'k:'); axis([0 T -25 25]);
xlabel('Time (s)'); ylabel('Response to Double Pulse (mm)')
