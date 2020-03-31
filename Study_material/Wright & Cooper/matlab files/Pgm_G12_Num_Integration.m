% Response of SDoF to single cycle of a square wave using SIMULINK
clear all; close all

% System parameters
f_nat = 2; period = 1 / f_nat; w_nat = 2 * pi * f_nat; zeta = 0.05;
mass = 1000; damp = 2 * zeta * mass * w_nat; 
stiff = w_nat^2 * mass; force = 1000;

% Simulation parameters
T = 8; dt = 0.01; t = [0:dt:T]; [dummy,nt] = size(t);

% Define excitation f(t)
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

% Run SIMULINK model for SDOF (works with column vectors in workspace 
% blocks Ein and Eout) using variable step length solver ODE45 with 
% outputs at same time intervals dt as for input – note that the SIMULINK file 
% needs to be in the same directory as the core MATLAB program
Ein = [t',f'];
sim('Model_G12_Num_Integration');

% Access stored data for response x - convert back to row vector for 
% consistency with f and t - convert to mm
x = 1000 * Eout';

% Plot response
figure(1); plot(t,x); axis([0 T -25 25]);
xlabel('Time (s)'); ylabel('Response to Single Cycle of a Square Wave');
