% Pgm_H2_Response_Gusts_Control 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up binary aeroelastic system and determines response to           % 
% control surface (chirp), gust (1 - cosine) and turbulence excitation   %
%                                                                        % 
% Parameters control_amp,  turb_amp, gust_amp_1_minus_cos define the     %
% amplitudes of the three possible inputs - set to zero to ignore        %
% tmax defines the maximum amount of data time sampled at dt             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = 100;    % AIRSPEED
% system parameters
s=7.5;         % semi span  
c=2;         % chord
a1 = 2*pi;    % lift curve slope
rho=1.225;    %air density      
m=100;          % unit mass / area of wing
kappa_freq = 5;  % flapping freq in Hz
theta_freq = 10; % pitch freq in Hz
xcm=0.5*c;   % position of centre of mass from nose 
xf=0.48*c;    % position of flexural axis from nose
Mthetadot=-1.2;   % unsteady aero damping term
e=xf/c - 0.25;  % eccentricity between flexural axis and aero centre

%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%
% control surface chirp
control_amp = 0;   % magnitude of control surface sweep input in degrees
c_burst = .333;   %  fraction of time length that is chirp signal  0 - 1
sweep_start = 1;  % chirp start freq in Hertz
sweep_end = 20;  % chirp end freq in Hertz

% "1 - Cosine" gust input 
gust_amp_1_minus_cos = 0;       % max velocity of "1 - cosine" gust   (m/s)
gust_t = .05;   %  fraction of total data time length that is gust    0 - 1 

% TURBULENCE INPUT - random signal between 0 Hz and turb_max_freq Hz
%   uniform amplitude across random frequency input
turb_amp = 10;       % max vertical velocity of turbulence   (m/s)
turb_t = 0.7;   %  fraction of total time length that is turbulence    0 - 1 
turb_max_freq = 20;    % max frequency of the turbulence (Hz) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  set up system matrices
a11=(m*s^3*c)/3 ;  % I kappa
a22= m*s*(c^3/3 - c*c*xf + xf*xf*c);   % I theta
a12 = m*s*s/2*(c*c/2 - c*xf);   %I kappa theta
a21 = a12;
k1=(kappa_freq*pi*2)^2*a11;     % k kappa
k2=(theta_freq*pi*2)^2*a22;     % k theta
A=[a11,a12;a21,a22];            % inertia matrix
E=[k1 0; 0 k2];                 % structural stiffness matrix
C=rho*V*[c*s^3*a1/6,0;-c^2*s^2*e*a1/4,-c^3*s*Mthetadot/8]; % aero damping
K=(rho*V^2*[0,c*s^2*a1/4;0,-c^2*s*e*a1/2])+[k1,0;0,k2] ;   % aero and structural stiffness

% gust vector
F_gust = rho*V*c*s*[s/4 c/2]';

% control surface vector
EE = .1;  % fraction of chord made up by control surface
ac = a1/pi*(acos(1-2*EE) + 2*sqrt(EE*(1-EE)));
bc = -a1/pi*(1-EE)*sqrt(EE*(1-EE));
F_control = rho*V^2*c*s*[-s*ac/4 c*bc/2]';

%  set up system matrices for SIMULINK 
MC=inv(A)*C;
MK=inv(A)*K;
MFG=inv(A)*F_gust;
MFC=inv(A)*F_control;

dt=0.001;   % sampling time
tmin=0;     % start time
tmax=10;    % end time
t=[0:dt:tmax]';       % Column vector

% CONTROL SURFACE INPUT SIGNAL - CHIRP SIGNAL
  control_amp = control_amp * pi / 180;   % radians
  t_end = tmax * c_burst;
  Scontrol = zeros(size(t));   % control input
  xt = sum(t < t_end);
if control_amp ~= 0
  for ii = 1:xt
    Scontrol(ii) = control_amp*sin(2*pi*(sweep_start + (sweep_end - sweep_start)*ii/(2*xt))*t(ii));
  end
else
    Scontrol(1:xt) = 0;
end

% GUST INPUT TERMS - "1-cosine"
  Sgust = zeros(size(t));
  %  1 - Cosine gust
  g_end = tmax * gust_t;
  gt = sum(t < g_end);
if gust_amp_1_minus_cos ~= 0
  for ii = 1:gt
    Sgust(ii) = gust_amp_1_minus_cos/2 * (1 - cos(2*pi*t(ii)/g_end));
  end
end

% TURBULENCE INPUT - unform random amplitude between 0 Hz and turb_max_freq Hz
  Sturb = zeros(size(t));
  t_end = tmax * turb_t;
  tpts = sum(t < t_end);  % number of turbulence time points required
  npts = max(size(t));
  if rem(max(tpts),2) ~= 0
    tpts =  tpts- 1;
  end
if turb_amp ~= 0
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
else
    Sturb(1:tpts)=zeros(tpts,1);
end
  for ii = 1:npts
    Sgust(ii) = Sgust(ii) + Sturb(ii); % "1 - cosine" plus turbulence inputs
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the system using SIMULINK
Egust=[t,Sgust]; % Gust Array composed of time and data columns
Econtrol=[t,Scontrol]; % Control Array composed of time and data columns

[tout] = sim('Model_H2_Response_Gusts_Control');

x1=EoutG1(:,1)*180/pi;   % kappa - flapping motion
x2=EoutG1(:,2)*180/pi;   % theta - pitching motion
x1dot=EoutG2(:,1)*180/pi;
x2dot=EoutG2(:,2)*180/pi;

figure(1)
plot(t,Scontrol,t,Sgust)
xlabel('Time (s)')
ylabel('Control Surface Angle(deg) and Gust Velocity(m/s)')
figure(2)
plot(t,x1,'r',t,x2,'b')
xlabel('Time (s)')
ylabel('Flap and Pitch Angles (deg/s)')
figure(3)
plot(t,x1dot,'r',t,x2dot,'b')
xlabel('Time (s)')
ylabel('Flap and Pitch Rates (deg/s)')
