% Add paths
added_paths

%Add sources & explanation per constant
constants

%Initial Position
rx = -4142343.75;
ry = 10*R_m;
R = [rx,ry,0];

%Initial Velocity
v = 7.1679e+03; %[m/s]
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R(1,:))^3*R(1,:);

%timesteps
dt_atmos = 0.5; %[s]
dt_kep_init = 1e-6; %[s]

%time to end simulation
tend = 3600 * 24 * 600; %[s]

%Control variables
control.a = 2.5*g_earth;
control.dalphadt = 1*pi/180;
control.dalpha = control.dalphadt*dt_atmos;
control.alpha_init = -25*pi/180; % rad
control.alpha_range = [-60 60]*pi/180;

% create atmosphere object
atm = marsatmosphere();
% create aerocoef object
aero_coef = aeroProperties('isotensoid');

