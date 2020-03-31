% Derivatives evaluated using Equivalent Air Speed and sea level air density

% Aerodynamic Derivatives (inertial axes) – Heave DoF
Z_0 = - 0.5 * rho0 * V0^2 * [- S_W * a_W + S_T * a_T * k_epsilon] * alpha_0;
Z_alpha = - 0.5 * rho0 * V0^2 * [S_W * a_W + S_T * a_T * (1 - k_epsilon)];
Z_q = - 0.5 * rho0 * V0 * S_T * a_T * l_T;
Z_eta = - 0.5 * rho0 * V0^2 * S_T * a_E;
Z_zdot = - 0.5 * rho0 * V0 * (S_W * a_W + S_T * a_T * (1 - k_epsilon));
Z_gW = - 0.5 * rho0 * V0 * S_W * a_W;
Z_gT = - 0.5 * rho0 * V0 * S_T * a_T * (1 - k_epsilon);

% Aerodynamic Derivatives (inertial axes) – Pitch DoF
M_0W = 0.5 * rho0 * V0^2 * S_W * c * C_M0 - 0.5 * rho0 * V0^2 * S_W * a_W * l_W * alpha_0;
M_0T = - 0.5 * rho0 * V0^2 * S_T * a_T * k_epsilon * l_T * alpha_0;
M_0 = M_0W + M_0T;
M_alpha = 0.5 * rho0 * V0^2 * [S_W * a_W * l_W - S_T * a_T * (1 - k_epsilon) * l_T];
M_q = - 0.5 * rho0 * V0 * S_T * a_T * l_T^2;
M_eta = - 0.5 * rho0 * V0^2 * S_T * a_E * l_T;
M_zdot = 0.5 * rho0 * V0 * (S_W * a_W * l_W - S_T * a_T * l_T * (1 - k_epsilon));
M_gW = 0.5 * rho0 * V0 * S_W * a_W * l_W;
M_gT = - 0.5 * rho0 * V0 * S_T * a_T * l_T * (1 - k_epsilon);

% Additional aerodynamic derivatives for wind axes – Rigid DoF (if required)
Z_w = - 0.5 * rho0 * V0 * (S_W * a_W + S_T * a_T * (1 - k_epsilon) + S_W * C_D);
M_w = 0.5 * rho0 * V0 * (S_W * a_W * l_W - S_T * a_T * l_T * (1 - k_epsilon));

% Aerodynamic Derivatives – Elastic DoF (if required)
Z_e = 0.5 * rho0 * V0^2 * (-S_W * a_W * J1 - S_T * a_T * gamma_eT);
Z_edot = - 0.5 * rho0 * V0 * S_T * a_T * kappa_eT;
M_e = 0.5 * rho0 * V0^2 * (S_W * a_W * l_W * J1 - S_T * a_T * l_T * gamma_eT);
M_edot = - 0.5 * rho0 * V0 * S_T * a_T * l_T * kappa_eT;
Q_0 = 0.5 * rho0 * V0^2 * (S_W * a_W * J2 - S_T * a_T * k_epsilon * kappa_eT) * alpha_0;
Q_alpha = 0.5 * rho0 * V0^2 * (-S_W * a_W * J2 - S_T * a_T * (1 - k_epsilon) * kappa_eT);
Q_q = - 0.5 * rho0 * V0 * S_T * a_T * l_T * kappa_eT;
Q_eta = -0.5 * rho0 * V0^2 * S_T * a_E * kappa_eT;
Q_e = 0.5 * rho0 * V0^2 * (- S_W * a_W * J3 - S_T * a_T * gamma_eT * kappa_eT);
Q_zdot = 0.5 * rho0 * V0 * (- S_W * a_W * J2 - S_T * a_T * (1 - k_epsilon) * kappa_eT);
Q_edot = - 0.5 * rho0 * V0 * S_T * a_T * kappa_eT^2;
Q_gW = - 0.5 * rho0 * V0 * S_W * a_W * J2;
Q_gT = - 0.5 * rho0 * V0 * S_T * a_T * kappa_eT;

% Additional aerodynamic derivative for wind axes – Elastic DoF
Q_w = 0.5 * rho0 * V0 * (- S_W * a_W * J2 - S_T * a_T * (1 - k_epsilon) * kappa_eT);
