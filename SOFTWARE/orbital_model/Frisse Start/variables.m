% Add paths
added_paths

%Add sources & explanation per constant
constants

%Initial Position
rx = -4148904.375000;
rx = -4655000;
ry = 10*R_m;
R = [rx,ry,0];

%Initial Velocity
v = 7.1679e+03; %[m/s]
v = 5.3757e3;
V = [0,-v,0];

%Initial acceleration
A = G*M_mars/norm(R(1,:))^3*R(1,:);

J = [0,0,0];

%timesteps
dt_atmos = 0.5; %[s]
dt_kep_init = 1e-6; %[s]

%time to end simulation
tend = 3600 * 24 *10; %[s]

%Control variables
control.a = 2.9*g_earth;
control.h = 60e3;
control.dalphadt = -1.0*pi/180;
control.dalpha = control.dalphadt*dt_atmos;
control.alpha_init = 15*pi/180; % rad
control.alpha_range = [-50 50]*pi/180;
control.control_up_lim = 100e3; % m
control.control_low_lim = 90e3; % m
control.Kp = 2*control.dalpha*10^(-6); % proportional gain
control.Ki = 5*control.dalpha*10^(-7); % integral gain
control.Kd = 2*control.dalpha*10^(-4); % differential gain
control.error = 0;
control.error_I = 0;


% create atmosphere object
atm = marsatmosphere();
% create aerocoef object
aero_coef = aeroProperties('irve');


