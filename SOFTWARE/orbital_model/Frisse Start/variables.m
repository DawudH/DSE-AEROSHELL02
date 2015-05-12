%Add sources & explanation per constant
constants

%Initial Position
ry = 2*R_m; %[m]
rx = -R_m-h_atm; %[m]
R = [rx,ry,0];

%Initial Velocity
v = 7000; %[m/s]
V = [0,-v,0];

%timesteps
dt_atmos = 0.5; %[s]
dt_kep_init = 1e-8; %[s]

%time to end simulation
tend = 3600 * 24 * 1; %[s]

%Control variables
control.CL_range = [-0.35 0.35];
control.CLCD = 0.25;
control.a = 2.9*g_earth;
control.CLa = 0.02;
control.dalpha = 0.2;
control.CL_init = -0.18;