% Data for Symmetric Aircraft
close all;  clear all

% Mass and Dimensions 
m = 10000; W = m * 9.81; m_F = 0.15 * m; m_W = 0.3 * m; m_C = 0.4 * m; m_T = 0.15 * m; 
S_W = 30.0; S_T = 7.5; s = 7.5; c = 2.0; s_tp = 3.0; c_tp = 1.25; l_W = 0.3*c; l_T = 3.5 * c; l_A = 0.125 * c; l_E = 0.125 * c; l_WM = l_W - l_A - l_E; l_M = 0.375 * c; l_F = (m_T * l_T - m_W * l_WM) / m_F; l_WT = l_W + l_T;
l_N = l_F; l_M = 0.375 * c;  mu = m_W / 2 / s;

% Moments of Inertia 
I_y_fuse = m_F * l_F^2 + m_T * l_T^2; I_y_W = m_W * (c / 3)^2;  I_y = I_y_fuse + I_y_W + m_W * l_WM^2;  l_y_W = sqrt(I_y_W / m_W); l_y = sqrt(I_y / m);

% Landing gear 
C_N = 3200; C_M = 19200; K_N = 80000; K_M = 240000; l_B = l_N + l_M;

% Basic aerodynamics 
a_W = 4.5; a_T = 3.2; a_E = 1.5; alpha_0 = - 0.03; C_M0 = - 0.03; C_D = 0.1; k_epsilon = 0.35; 
