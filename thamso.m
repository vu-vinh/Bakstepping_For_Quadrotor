%% Process
% thamso -> SMC -> thamso(reset parameters) -> BS -> dothi

%% init value
clear; clearvars;

global m g Ixx Iyy Izz Ir k b l Ax Ay Az M
global d w i p_d p_d_dot p_d_2dot a_d a_d_dot a_d_2dot
global T tau_phi tau_theta tau_psi ux uy uz Ts

m = 1.18;         % Quadrotor mass (kg)
g = 9.81;         % Gravity acceleration (m/s^2)

Ixx = 0.01;       % Moments of inertia (kg*m^2)
Iyy = 0.01;
Izz = 0.0148;
Ir = 2.83e-5;     % Rotor moment of inertia

k = 2.98e-5;      % Lift constant
b = 3.23e-7;      % Drag constant
l = 0.25;         % Arm length (m)

Ax = 5.56e-4;     % Drag coefficients
Ay = 5.56e-4;
Az = 6.35e-4;

M = [1/(4*k)  0          -1/(2*k*l) -1/(4*b);
     1/(4*k) -1/(2*k*l)   0          1/(4*b);
     1/(4*k)  0           1/(2*k*l) -1/(4*b);
     1/(4*k)  1/(2*k*l)   0          1/(4*b)];

% Step
Ts = 0.01;
Tf = 45;
N = Tf/Ts + 1;
t = 0:Ts:Tf;

%% Initial state
% State order: [x y z phi theta psi vx vy vz phi_dot theta_dot psi_dot]'
t0 = 0;
x0 = zeros(12,1);
% x0(3) = 2;

ux = 0;
uy = 0;
uz = 0;

xout = zeros(12,N);
xout(:,1) = x0;
u = zeros(4,N);  
w = zeros(4,N); 


%% Scenario 1: elipse
xd = 0.5*cos(pi*t/20);
yd = 0.5*sin(pi*t/20);
zd = 2 - 0.5*cos(pi*t/20);

xd_dot = -(0.5*pi/20)*sin(pi*t/20);
yd_dot = (0.5*pi/20)*cos(pi*t/20);
zd_dot = (0.5*pi/20)*sin(pi*t/20);

xd_2dot = zeros(1,N); 
yd_2dot = zeros(1,N); 
zd_2dot = zeros(1,N);

% w0 = pi/20;
% xd = 2*sin(w0*t).*cos(w0*t/2);
% yd = 2*sin(w0*t).*sin(w0*t/2);
% zd = 2 - sin(w0*t).*(cos(w0*t/2) + sin(w0*t/2));
% 
% % Dao ham bac 1
% xd_dot = (3*w0/2)*cos(3*w0*t/2) + (w0/2)*cos(w0*t/2);
% yd_dot = (3*w0/2)*sin(3*w0*t/2) - (w0/2)*sin(w0*t/2);
% 
% zd_dot = -w0*cos(w0*t).*(cos(w0*t/2) + sin(w0*t/2)) - (w0/2)*sin(w0*t).*(cos(w0*t/2) - sin(w0*t/2));
% 
% % Dao ham bac 2
% xd_2dot = -(9*w0^2/4)*sin(3*w0*t/2) - (w0^2/4)*sin(w0*t/2);
% 
% yd_2dot = (9*w0^2/4)*cos(3*w0*t/2) - (w0^2/4)*cos(w0*t/2);
% 
% zd_2dot = (5*w0^2/4)*sin(w0*t).*(cos(w0*t/2) + sin(w0*t/2)) - w0^2*cos(w0*t).*(cos(w0*t/2) - sin(w0*t/2));

% Disturbances
d = zeros(6,N);

d(1,:) = 0.3*sin(t);       % dx
d(2,:) = 0.3*sin(t);       % dy
d(3,:) = 0.1*cos(0.5*t);   % dz
d(4,:) = 0.2*sin(2*t);     % dphi
d(5,:) = 0.2*sin(2*t);     % dtheta
d(6,:) = 0.1*sin(t);       % dpsi

p_d = [xd; yd; zd];
p_d_dot = [xd_dot; yd_dot; zd_dot];
p_d_2dot = [xd_2dot; yd_2dot; zd_2dot];

a_d = zeros(3,N);
a_d_dot = zeros(3,N);
a_d_2dot = zeros(3,N);
