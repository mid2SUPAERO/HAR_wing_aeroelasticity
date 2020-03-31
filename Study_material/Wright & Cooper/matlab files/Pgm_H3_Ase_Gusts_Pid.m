% Pgm_H3_Ase_Gusts_Pid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up binary aeroelastic system    %%
% Determines the response at one speed from gust(1 - cosine)
% and / or turbulence input
% Applies PID control with Kv and Kd gains for velocity and displacement 
% terms via a control surface
% Assumes that the transducer is at the wing tip leading edge.
% Calls simulink routine   Model_H3_Ase_Gusts_Pid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%
% System Parameters
s=7.5;         % semi span
c=2;         % chord
a = 2*pi;    % lift curve slope
m=100;          % unit mass / area of wing
rho=1.225;    %air density      
kappa_freq = 5;  % bending (flapping) freq in Hz
theta_freq = 10; % torsion (pitch) freq in Hz
xcm=0.5*c;   % position of centre of mass from nose 
xf=0.48*c;    % position of flexural axis from nose
e=xf/c - 0.25;  % eccentricity between flexural axis and aero centre (1/4 chord)
Mthetadot=-1.2;
V = 100;    % velocity system is considered

%   GAINS
kv = -.01;   % gain for velocity terms
kd = 0;   % gain for displacement terms

%  1 - COSINE GUST PARAMETERS
gust_amp_1_minus_cos = 0;       % max velocity of "1 - cosine" gust   (m/s)
gust_t = .1;   %  fraction of total time length that is gust    0 - 1 
% TURBULENCE INPUT PARAMETERS
turb_amp = 10;       % max vertical velocity of turbulence   (m/s)
turb_t = .8;   %  fraction of total time length that is turbulence    0 - 1 
turb_max_freq = 20;    % max frequency of the turbulence (Hz) 


% system parameters     
a11=(m*s^3*c)/3;  % I kappa
a22= m*s*(c^3/3 - c*c*xf + xf*xf*c); % I theta
a12 = m*s*s/2*(c*c/2 - c*xf);  %I kappa theta
a21 = a12;

k1=(kappa_freq*pi*2)^2*a11;  %k kappa
k2=(theta_freq*pi*2)^2*a22;  %k theta

A=[a11,a12;a21,a22];
E=[k1 0; 0 k2];

EE = .1;  % fraction of chord made up by control surface
ac = a/pi*(acos(1-2*EE) + 2*sqrt(EE*(1-EE)));
bc = -a/pi*(1-EE)*sqrt(EE*(1-EE));

% control surface terms
F_control = rho*V^2*c*s*[-s*ac/4 c*bc/2]';
g1 = F_control(1);
g2 = F_control(2);

F = kv*[g1*s -g1*c/2;  g2*s -g2*c/2];  % control damping matrix
H = kd*[g1*s -g1*c/2;  g2*s -g2*c/2];  % control damping matrix

% include the control terms into the damping and stiffness matrices
C = rho*V*[c*s^3*a/6,0;-c^2*s^2*e*a/4,-c^3*s*Mthetadot/8] - F;
K=(rho*V^2*[0,c*s^2*a/4;0,-c^2*s*e*a/2])+[k1,0;0,k2] - H ;

% gust terms (rhs of aeroelastic equations)
F_gust = rho*V*c*s*[s/4 c/2]';

%  set up system matrices
MC=inv(A)*C;
MK=inv(A)*K;
MFG=inv(A)*F_gust;

dt=0.001;
tmin=0;
tmax=10;
t=[0:dt:tmax]';       % Column vector
npts = max(size(t));

% GUST INPUT TERMS
Sgust = zeros(size(t));
%  1 - Cosine gust
g_end = tmax * gust_t;
gt = sum(t < g_end);
if gust_amp_1_minus_cos ~= 0
  for ii = 1:gt
    Sgust(ii) = gust_amp_1_minus_cos/2 * (1 - cos(2*pi*t(ii)/g_end));
  end
end

% TURBULENCE INPUT
Sturb = zeros(size(t));
if turb_amp ~= 0
  t_end = tmax * turb_t;
  tpts = sum(t < t_end);  % number of turbulence time points required


  if rem(max(tpts),2) ~= 0
    tpts = tpts - 1;
  end
  td2 = tpts / 2;
  td2p1 = tpts/2 + 1;
  df = 1/(tpts*dt);
  fpts = fix(turb_max_freq / df) + 1; %  number of frequency points that are required to form turbulence series

  for ii = 1:fpts              % define real and imag parts of freq domain terms magnitude of unity, random phase
    a(ii) = 2 * rand - 1;    % real part  - 1 < a < 1
    b(ii) = sqrt(1 - a(ii)*a(ii)) * (2*round(rand) - 1);   % imag part  
  end
    
  tf =  (a + j*b);   %  determine complex frequency representation with correct frequency characteristics
  tf(fpts+1 : td2p1) = 0;
  tf(td2p1+1 : tpts) = conj(tf(td2:-1:2));
  Sturb(1:tpts) = turb_amp * real(ifft(tf));
  kscale = turb_amp / max(Sturb);
  Sturb = kscale*Sturb;
end

for ii = 1:npts
    Sgust(ii) = Sgust(ii) + Sturb(ii);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second Order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Egust=[t,Sgust]; % Gust Array composed of time and data columns

[tout] = sim('Model_H3_Ase_Gusts_Pid');

kappa=EoutG1(:,1)*180/pi;   % kappa - flapping motion
theta=EoutG1(:,2)*180/pi;   % theta - pitching motion
kappadot=EoutG2(:,1)*180/pi;
thetadot=EoutG2(:,2)*180/pi;

% determine the control response
beta = kv*(kappadot * s - thetadot*c/2) + kd*(kappa * s - theta*c/2);  % control surface (in DEGREES)

% determine the measurement position response
z = kappa * s - theta * c/2;

figure(1)
plot(t,Sgust)
xlabel('Time (s)')
ylabel('Gust Velocity (m/s)')

figure(2)
subplot(211)
plot(t,kappa,'k')
xlabel('Time (s)')
ylabel('Kappa Response (deg)')
grid
subplot(212)
plot(t,theta,'k')
xlabel('Time (s)')
ylabel('Theta Response (deg)')
grid

figure(3)
plot(t,kappadot,'r',t,thetadot,'b')
xlabel('Time (s)')
ylabel('System Response (deg/s)')


figure(4)
plot(t,beta)
xlabel('Time (s)')
ylabel('Control Angle (deg)')
grid

figure(5)
plot(t,z)
xlabel('Time (s)')
ylabel('Wing Tip LE Displacement (m)')
grid


