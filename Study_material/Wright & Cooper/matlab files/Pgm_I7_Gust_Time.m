% Gust response of a rigid aircraft in the time domain – requires aircraft data, flight case and derivative codes
close all;
%Gust profile                                                                  ***************
%Air space   ---------------------------------------------------------------------------------------
%Aircraft      W+++++++T                                                                 W+++++++T

% Set-up time array and length of simulation
tmin = 0; tmax = 8.0; 
dt = 0.005;  		% Time increment - will need to be smaller for the elastic aircraft (~ 0.002s)
t = [tmin: dt: tmax]';  [Nsim, dummy] = size(t);   

% Set-up distance array and size for simulation (important to use TAS not EAS)
x_g = V * t; 			    				% Distance array
dx = V * dt;			    				% Distance increment
Nb = round(l_WT / dx+1);                     				% No of points between wing / tailplane
N = Nsim + Nb - 1;                                				% No of points for total air space

% Gust velocity profile (1 – cosine)
delta_wgt = 5.0;                                     				% Max gust velocity TAS
L_g = 250.0;			     				% Gust length
Nd = round(L_g / dx + 1);                      				% No of points for gust
Nbd = Nb + Nd - 1;                                				% No of points for a/c and gust
wg = zeros(N, 1); wg_W = zeros(N, 1); wg_T = zeros(N, 1);
wg(Nb:Nbd) = (delta_wgt / 2) * (1 - cos(2 * pi * x_g(1:Nd) / L_g));    	% Gust velocity array at centre of mass
wg_W(1:Nsim) = wg(Nb:N); wg_T(1:Nsim) = wg(1:Nsim);                	% Gust velocity array at wing and tailplane

% Equations of motion for rigid aircraft (elastic is similar but with an additional row / column for elastic mode
M = [m  0; 0  I_y]; CC = - [Z_zdot  Z_q; M_zdot  M_q]; KK = - [0  Z_alpha; 0  M_alpha];
FW = [Z_gW;  M_gW]; FT = [Z_gT;  M_gT];
MI = inv(M); MC = MI * CC; MK = MI * KK; MFW = MI * FW; MFT = MI * FT;

% Simulation to find displacements and velocities of response in inertial axes
inW = [t, wg_W(1:Nsim)]; inT = [t, wg_T(1:Nsim)];                             	% Input arrays at wing and tail[lane
sim('Model_I7_Gust_Time')                                              % Solve equations via SIMULINK model

% Displacements, velocities and accelerations at CoM, nose and tail (all Nsim x 1 arrays)
z_C = out_disp(:,1);  zdot_C = out_rate(:,1);  zddot_C = out_acc(:,1);  
theta = out_disp(:,2);  theta_dot = out_rate(:,2);  theta_ddot = out_acc(:,2);
z_F = z_C - l_F * theta;  zddot_F = zddot_C - l_F * theta_ddot;  
z_T = z_C + l_T * theta;  zddot_T = zddot_C + l_T * theta_ddot;

% Plot responses
subplot(211); plot(t, theta * 180 / pi, 'k-', t, theta_dot * 180 / pi, 'k:')
title('Pitch Response for Rigid Aircraft in Gust')
xlabel('Time (s)'); ylabel('Pitch (deg) and Pitch Rate (deg/s)'); legend('Pitch Angle', 'Pitch Rate')
subplot(212); plot(t, zddot_C / 9.81, 'k-', t, zddot_T / 9.81, 'k:')
title('Nose / CoM / Tail Heave Acceleration for Rigid Aircraft in Gust')
xlabel('Time (s)'); ylabel('Heave Acceleration (g)'); legend('CoM', 'Tailplane')
