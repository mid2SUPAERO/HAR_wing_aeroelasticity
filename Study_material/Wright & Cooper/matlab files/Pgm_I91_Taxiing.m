% Rigid aircraft taxiing over a 1 – cosine dip – requires aircraft data, flight case but not derivative codes
close all;
% Time and distance data
tmin = 0; tmax = 10.0; dt = 0.01;  % Time increment - will need to be smaller for the elastic aircraft (~ 0.002s)
t = [tmin: dt: tmax]';  [Nsim, dummy]=size(t);
x_r = V * t; dx = V * dt; Nb = round(l_B / dx + 1); N = Nsim + Nb - 1;
% Runway profile – define (1 – cosine) dip – note similarity to gust except for the need for hdot terms
 
%Dip                                                                       ***************
%Runway   ----------------------------------------------------------------------------------
%Aircraft    N+++++++M                                                            N+++++++M
 
delta_h_r = 0.03;  L_r = 60;  Ndip = round(L_r / dx + 1); Nbd = Nb + Ndip - 1;
h = zeros(N, 1);  h_N = zeros(N, 1);  h_M = zeros(N, 1);  
hdot = zeros(N, 1);  h_Ndot = zeros(N, 1); h_Mdot = zeros(N, 1);  
h(Nb:Nbd) = (delta_h_r / 2) * (1 - cos(2 * pi * x_r(1:Ndip) / L_r));
hdot(Nb:Nbd) = V * pi * delta_h_r / L_r * sin(2 * pi * x_r(1:Ndip) / L_r);
 
% Runway profile defined at nose and main gear positions
h_N(1:Nsim) = h(Nb:N);  h_M(1:Nsim) = h(1:Nsim);  
h_Ndot(1:Nsim) = hdot(Nb:N);  h_Mdot(1:Nsim) = hdot(1:Nsim);
 
% Equations of motion in second order form
M = [m  0; 0  I_y];
CC = [C_N + C_M   -l_N * C_N + l_M * C_M;   -l_N * C_N + l_M * C_M   l_N^2 * C_N + l_M^2 * C_M];
K = [K_N + K_M   -l_N * K_N + l_M * K_M;    -l_N * K_N + l_M * K_M    l_N^2 * K_N + l_M^2 * K_M];
DC = [C_N  C_M;   -l_N * C_N   l_M*C_M];  DK = [K_N   K_M;  -l_N * K_N   l_M * K_M];
MC = inv(M) * CC;  MK = inv(M) * K;  MDC = inv(M) * DC;  MDK = inv(M) * DK;
 
% State space equations in first order form
Null22 = zeros(2, 2);  I = eye(2);  C = eye(4); D = zeros(4, 4);
A = [Null22  I;  -MK  -MC];  B = [Null22  Null22;  MDK  MDC];
 
% State space input vector (Nsim x 4)
u = [h_N(1:Nsim)  h_M(1:Nsim)  h_Ndot(1:Nsim)  h_Mdot(1:Nsim)];
 
% Simulation in state space to find displacements and velocities
in = [t, u];  sim('Model_I91_Taxiing')
 
% Displacements at CoM, nose and main gears (all Nsim x 1)
% Dimension of 'out' is Nsim x 4 (zc, theta, zcdot, thetadot)
z_C = out(:,1);  z_N = z_C - out(:,2) * l_N;  z_M = z_C + out(:,2) * l_M;
 
% Calculate accelerations from state space equations
% Dimension of 'outdot' is Nsim x 4 (zcdot, thetadot, zcdddot, thetaddot)
outdot = (A * out' + B * u')'; zddot_C = outdot(:,3);  theta_ddot = outdot(:,4); 
zddot_N = zddot_C - theta_ddot * l_N;  zddot_M = zddot_C + theta_ddot * l_M;
 
% Plot responses at nose and main gear
subplot(211); plot(t, 1000 * z_N, 'k-', t, 1000 * z_M, 'k:')
title('Nose / Main Gear Heave for Taxiing Rigid Aircraft')
xlabel('Time (s)'); ylabel('Heave Response (mm)'); legend('Nose Gear', 'Main Gear')

subplot(212); plot(t, zddot_N / 9.81, 'k-', t, zddot_M / 9.81,'k:')
title('Nose / Main Gear Heave Acceleration for Taxiing Rigid Aircraft')
xlabel('Time (s)'); ylabel('Heave Acceleration (g)'); legend('Nose Gear', 'Main Gear')
