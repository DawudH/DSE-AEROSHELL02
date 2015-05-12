% Add paths
added_pathes

%Add sources & explanation per constant
constants

%Initial Position
ry = SOI; %[m]
rx = -4190000; %[m]
R = [rx,ry,0];

%Initial Velocity
v = 7000; %[m/s]
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R)^3*R;

%timesteps
dt_atmos = 0.5; %[s]
dt_kep_init = 1e-6; %[s]

%time to end simulation
tend = 3600 * 24 * 1; %[s]

%Control variables
control.CL_range = [-0.35 0.35];
control.CLCD = 0.25;
control.a = 2.9*g_earth;
control.CLa = 0.02;
control.dalpha = 0.2;
control.CL_init = -0.18;

% create atmosphere object
atm = marsatmosphere();