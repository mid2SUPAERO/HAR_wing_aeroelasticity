% Dynamic response of a symmetric rigid aircraft to an elevator input – requires aircraft data, flight case and derivative codes
close all;
% Simulation data
tmin = 0;  tmax = 5.0;  dt = 0.005;       % Time increment - will need to be smaller for the elastic aircraft (~ 0.002s)
t = [tmin: dt: tmax]';  [N, dummy] = size(t);  eta = zeros(N, 1);

% On / off pulse input to elevator
tpulse = 1.0;  npulse = tpulse / dt + 1;  eta_in = 2.0;  eta_in = eta_in * pi / 180;
eta(1:npulse) = eta_in * ones(npulse, 1);

% Aircraft initial condition
U_e = V0;  W_e = 0;

% Flight mechanics linearised equations of motion
M = [m   0; 0   I_y];  C = [-Z_w  -(m * U_e + Z_q);  -M_w  -M_q];  F = [Z_eta;  M_eta];
MC = inv(M) * C;  MF = inv(M) * F;

% Simulation of response in body fixed (wind) axes to yield w,q and then theta
in = [t, eta];
sim('Model_I6_Dynamic_Manoeuvres');

% Output variables
w = out_rate(:,1);						% Downwards velocity
q = out_rate(:,2);  						% Pitch rate
wdot = out_acc(:,1); 					% Downwards acceleration
theta = out_disp(:,2);            				% Integral of q (i.e. theta)
alpha = w / V0;                 					% Incidence
gamma = theta - alpha;              				% Flight path angle - perturbation
az = wdot - q * U_e;                  				% Normal acceleration at CoM

% Response plots – body fixed axes
figure(1);  subplot(311);  plot(t, q * 180 / pi, 'k-', t, alpha * 180 / pi, 'k:')
title('Pitch Rate and Incidence Response for Elevator Input')
xlabel('Time (s)');  ylabel('Pitch Rate (deg/s), Incidence (deg)');  legend('Pitch Rate', 'Incidence')
subplot(312);  plot(t, az / 9.81, 'k-')
title('Normal Acceleration for Elevator Input');  xlabel('Time (s)');  ylabel('Normal Acceleration (g)')
subplot(313);  plot(t, theta * 180 / pi, 'k-', t, gamma * 180 / pi, 'k:')
title('Pitch and Flight Path Angles for Elevator Input')
xlabel('Time (s)');  ylabel('Pitch and Flight Path Angles (deg)');  legend('Pitch Angle', 'Flight Path Angle')

% Transform velocities related to CoM from body fixed to earth fixed axes
% Centre of Mass
U_E = U_e + W_e * theta;
W_E = -U_e * theta + W_e + w;
% Tailplane
U_E_tp = U_e + W_e * theta + w.* theta;
W_E_tp = -U_e * theta + W_e + w + l_T * q;

% Response plots – velocities in earth axes
figure(2);  plot(t, w, 'k-', t, W_E, 'k:')
title('Vertical Velocity Response (relative to wind / earth axes) for Elevator Input')
xlabel('Time (s)');  ylabel('Velocity (m/s)')
legend('Vertical velocity relative to wind axes w', 'Vertical velocity relative to earth axes WE')
figure(3);  plot(t, W_E, t, W_E_tp, 'k-')
title('Vertical Velocity Response (relative to earth axes) for Elevator Input')
xlabel('Time (s)');  ylabel('Velocity WE (m/s)');  legend('Centre of mass', 'Tailplane')

% Integrate to yield CoM position coordinates from velocities in earth axes
sim('Model_I6_Earth_Axes');

% Response plots – CoM position in earth axes
figure(4);  plot(X_E, -Z_E, 'k-')
title('CoM Flight Profile in Earth Axes following Elevator Input')
xlabel('Horizontal Displacement  XE (m)');  ylabel('Vertical Displacement  ZE (m)')
