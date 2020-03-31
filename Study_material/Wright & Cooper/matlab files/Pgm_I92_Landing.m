% Rigid aircraft landing - non-linear gas spring / tyre model - requires aircraft data, but not flight case / derivatives
close all;
mh = m / 2;  					% Factor total mass for half aircraft
Wh = mh * 9.81;
Lh = Wh;					% Lift balances weight
W_e = 3.0;  					% Vertical speed on landing (TAS)

% Single main landing gear parameters leading to non-linear behaviour
Pinf = 25e5;  PS = 1e7;  PC = 3e7;  PA = 1e5;  zS = 0.4;  Area = Wh / PS;
Vratio = PC / Pinf;  Vinf = Vratio * Area * zS / (Vratio - 1);
zinf = Vratio / (Vratio - 1) * zS;  ns = 1;  nd = 1.35;

% Generate force ~ displacement variation for non-linear shock absorber stiffness look-up table
dz = 0.005;  z = [0: dz: zS];  [dummy, nz] = size(z);
for j=1:nz
    Pd(j) = Pinf / (1 - z(j) / zinf)^nd;
    Fd(j) = (Pd(j) - PA) * Area;
end

% Generate force ~ velocity variation for non-linear shock absorber damping look-up table
Ccomp = 8000;  Crecoil = 120000;
zdot = [-1.5: 0.005: 3];  [dummy, nzdot] = size(zdot);
for j=1:nzdot
    if zdot(j) >= 0
        Fdamp(j) = Ccomp * zdot(j)^2;
    else
        Fdamp(j) = - Crecoil * zdot(j)^2;
    end
end

% Tyre data for half aircraft
mt = 100;  kt = 1000e3; ct = 0;

% Run simulation for aircraft half mass on non-landing gear with tyre mass / spring representation
tmin = 0;  tmax = 0.5;  dt = 0.0002;  t = [tmin: dt: tmax]';
sim('Model_I92_Landing')

% Plot output
figure(1);  plot(t, -out_accnl / 9.81,'k-');  title('Main Gear Deceleration for Rigid Aircraft Landing - Non-linear')
xlabel('Time (s)');  ylabel('Main Gear Deceleration (g)')

figure(2);  plot(t, out_dispnl, 'k-', t, (out_dispnl - out_dispnlt), 'k:', t, out_dispnlt, 'k--')
title('Main Gear / Shock Absorber / Tyre Displacement Response of Rigid Aircraft Landing - Non-linear')
xlabel('Time (s)');  ylabel('Displacement (m)');  legend('Main Gear', 'Shock Absorber', 'Unsprung Mass')

Fgrnd = ct * out_ratenlt + kt * out_dispnlt;		% Ground reaction force
figure(3);  plot(t, Fgrnd / Wh, 'k-');  title('Normalised Ground Load - Rigid Aircraft Landing - Non-linear')
xlabel('Time (s)');  ylabel('Normalised Ground Load')
