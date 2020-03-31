% Data for Flight Case – specify EAS (and relative density if not at sea level)
V0 = 150;  			% EAS 
rho0 = 1.225;			% Sea level density
root_sigma  = 0.8;        		% Relative density at altitude
V = V0 / root_sigma;                 	% Convert to TAS for gusts / ground manoeuvres
