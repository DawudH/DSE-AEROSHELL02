% Add paths
added_paths

%Add sources & explanation per constant
constants

%initial angles
theta0 = 135/180*pi;%3.016549;%150/180*pi;
%d=12
gamma = 21.88;
% %d=6
% gamma = 22.2;
%d=18
%gamma = 21.74;

%Initial Position
rx = -4148904.375000;
rx = -4655000;
ry = 10*R_m;
R = [rx,ry,0];
r = h_atm+R_m;

%Initial Velocity
v = 7.1679e+03; %[m/s] at 10Rm for 7km/s at SOI
v = 5.3757e3; %[m/s] at 10Rm for 7km/s at atmos
v = 7000; %[m/s] at atmos for 7km/s at atmos
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R(1,:))^3*R(1,:);

J = [0,0,0];

%timesteps
dt_atmos = 0.1; %[s]
dt_kep_init = 1e-6; %[s]

%time to end simulation
tend = 3600 * 24 *10; %[s]

%Control variables
control.a = 2.9*g_earth;
control.h = 60e3;
control.dalphadt = -0.1*pi/180;
control.dalpha = control.dalphadt*dt_atmos;
control.alpha_init = 10*pi/180; % rad
control.alpha_range = [-50 50]*pi/180;
control.control_up_lim = 100e3; % m
control.control_low_lim = 90e3; % m
control.Kp = 2*control.dalpha*10^(-6); % proportional gain
control.Ki = 5*control.dalpha*10^(-7); % integral gain
control.Kd = 2*control.dalpha*10^(-4); % differential gain
control.error = 0;
control.error_I = 0;

Crho = 1.0;

% create atmosphere object
atm = marsatmosphere();
% create aerocoef object
aero_coef = aeroProperties('irve');


