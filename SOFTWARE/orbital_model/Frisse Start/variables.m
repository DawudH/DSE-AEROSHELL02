% Add paths
added_paths

%Add sources & explanation per constant
constants

%Initial Position
rx = -4146953.125000;
ry = 10*R_m;
R = [rx,ry,0];

%Initial Velocity
v = 7.1679e+03; %[m/s]
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R(1,:))^3*R(1,:);

J = [0,0,0];

%timesteps
dt_atmos = 0.5; %[s]
dt_kep_init = 1e-6; %[s]

%time to end simulation
tend = 3600 * 24 * 5; %[s]

%Control variables
control.a = 2.5*g_earth;
control.dalphadt = 1*pi/180;
control.dalpha = control.dalphadt*dt_atmos;
control.alpha_init = 10*pi/180; % rad
control.alpha_range = [-25 25]*pi/180;
control.Kp = 1; % proportional gain

% create atmosphere object
atm = marsatmosphere();
% create aerocoef object
aero_coef = aeroProperties('apollo');

