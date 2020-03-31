% Pgm_H1_Calcs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up the aeroelastic matrices for binary aeroelastic model          %
% performs eigenvalue solution at desired speeds and determines the      %
% frequencies and dampings ratios                                        %
% plots v_omega and v_g plots trends                                     %
%                                                                        %
% vary the speed range via velstart and velend parameters                %
% kappa_freq and theta_freq are flap and pitch nat freqs (Hz)            %
% damping_Y_N defines whether damping (aero and structural) is included  %
% z1 and z2 are structural damping ratios for both modes (if required)   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize variables
clear
clf
s=7.5;          % semi span
c=2;            % chord
m=100;          % unit mass / area of wing
kappa_freq = 5; % flapping freq in Hz
theta_freq = 10;% torsion (pitch) freq in Hz
xcm=0.5*c;     % position of centre of mass from nose 
xf=0.48*c;      % position of flexural axis from nose
e=xf/c - 0.25;  % eccentricity between flexural axis and aero centre (1/4 chord)

velstart = 1;   % lowest velocity to consider
velend = 240;   % maximum velocity
velinc=.1;      % velocity increment
a=2*pi;         % 2D lift curve slope
rho=1.225;      %air density      
Mthetadot=-1.2;   % unsteady damping term
M = (m*c^2 - 2*m*c*xcm)/(2*xcm);   % leading edge mass term to allow mass axis to move

damping_Y_N = 1;  % =1 if damping  included   =0 if not included
if damping_Y_N == 1
    % structural proportional damping inclusion    C = alpha * M + beta * K
    % two freqs and damps must be defined
% set dampings to zero for no structural damping
    z1 = 0.0;               % critical damping at first frequency
    z2 = 0.0;               % critical damping at second frequency
    w1 = 2*2*pi;            % first frequency
    w2 = 14*2*pi;           % second frequency
    alpha = 2*w1*w2*(-z2*w1 + z1*w2)/ (w1*w1*w2*w2);
    beta = 2*(z2*w2-z1*w1) / (w2*w2 - w1*w1);
end

% set up system matrices =======================================
% inertia matrix
a11=(m*s^3*c)/3 + M*s^3/3;  % I kappa
a22= m*s*(c^3/3 - c*c*xf + xf*xf*c) + M*(xf^2*s); % I theta
a12 = m*s*s/2*(c*c/2 - c*xf) - M*xf*s^2/2;  %I kappa theta
a21 = a12;
A=[a11,a12;a21,a22];
% structural stiffness matrix
k1=(kappa_freq*pi*2)^2*a11; %k kappa
k2=(theta_freq*pi*2)^2*a22; %k theta
E=[k1 0; 0 k2];

icount = 0;
for V=velstart:velinc:velend       % loop for different velocities ===========================
    icount = icount +1;
    if damping_Y_N == 0;    % damping matrices
        C=[0,0;0,0];        % =0 if damping not included
    else                    % =1 if damping included
        C=rho*V*[c*s^3*a/6,0;-c^2*s^2*e*a/4,-c^3*s*Mthetadot/8] + alpha*A + beta*E; % aero and structural damping
    end
    K=(rho*V^2*[0,c*s^2*a/4;0,-c^2*s*e*a/2])+[k1,0;0,k2];  % aero and structural stiffness

    Mat=[[0,0;0,0],eye(2);-A\K,-A\C];   % set up 1st order eigenvalue solution matrix
    lambda = eig(Mat);           % eigenvalue solution

    %natural frequencies and damping ratios
    for jj = 1:4
        im(jj)=imag(lambda(jj));
        re(jj)=real(lambda(jj));
        freq(jj,icount)=sqrt(re(jj)^2+im(jj)^2);
        damp(jj,icount)=-100*re(jj)/freq(jj,icount); 
        freq(jj,icount) = freq(jj,icount)/(2*pi);    % convert frequency to hertz
    end
    Vel(icount)=V;
end

% Plot frequencies and damping ratios vs speed
figure(1)
subplot(2,1,1); plot(Vel,freq,'k');
vaxis = axis; xlim = ([0 vaxis(2)]);
xlabel ('Air Speed (m/s) ' ); ylabel ('Freq (Hz)'); grid

subplot(2,1,2);
plot(Vel,damp,'k')
xlim = ([0 vaxis(2)]); axis([xlim ylim]);
xlabel ('Air Speed (m/s) ' ); ylabel ('Damping Ratio (%)'); grid






