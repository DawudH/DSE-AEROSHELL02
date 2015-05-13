% Add paths
added_paths

%Add sources & explanation per constant
constants

%Initial Position
step_rx = -10000;
rx_b = -4100000;
rx_e = -4500000;
rx = (rx_b:step_rx:rx_e)'; %[m]
n = length(rx);
ry = SOI*ones(n,1); %[m]
R = [rx,ry,zeros(n,1)];
% rx = -4190000;
% ry = SOI;
% R = [rx,ry,0];
%Initial Velocity
v = 7000; %[m/s]
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R(1,:))^3*R(1,:);

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