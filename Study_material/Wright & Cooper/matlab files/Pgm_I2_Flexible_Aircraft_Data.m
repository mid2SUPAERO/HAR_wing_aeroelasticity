% Additional Data for Flexible Aircraft (if required) - Example of Fuselage Bending dominant
kappa_e0 = 1; gamma_e0 = 0; A = 0; B = 0;
    
% Normalisation to 1 at wing tip trailing edge
kappa_tip = kappa_e0 * (1 + A) + gamma_e0 * (1 + B) * (c - c/4 - l_A);
kappa_e0 = kappa_e0 / kappa_tip;

% Solution of equations to yield dominant Fuselage Bending mode shape
X = [m_F  m_T;  -m_F * l_F  m_T * l_T]; Y = [- (m_W + m_C) * kappa_e0;  m_W * l_WM * kappa_e0]; Z = X\Y;
kappa_eC = 1; kappa_eF = Z(1); kappa_eT = Z(2); gamma_eT = 2 * (kappa_eT - kappa_eC) / l_T - gamma_e0;

% Modal Mass, Natural Frequency and Modal Stiffness for Fuselage Bending dominant
m_e = m_F * kappa_eF^2 + m_W * kappa_e0^2 + m_C * kappa_eC^2 + m_T * kappa_eT^2;
f_e = 4.0;  omega_e = 2 * pi * f_e;  k_e = omega_e^2 * m_e;

% ‘J’ Integrals for aerodynamic derivatives
J1 = gamma_e0 * (1 + B / 2);
J2 = kappa_e0 * (1 + A / 3) - l_A * gamma_e0 * (1 + B / 2);
J3 = gamma_e0 * kappa_e0 * (1 + A / 3 + B / 2 + A * B / 4) - l_A * gamma_e0^2 * (1 + B + B^2 / 3);
