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
control.dalphadt = -0.2*pi/180;
control.dalpha = control.dalphadt*dt_atmos;
control.alpha_init = 20*pi/180; % rad
control.alpha_range = [-50 50]*pi/180;
control.control_up_lim = 2; % g_earth
control.control_low_lim = 0.5; % g_earth
control.Kp = 1*0.005*control.dalpha; % proportional gain
control.Ki = 1* 0.0001*control.dalpha; % integral gain
control.Kd = 1* 2*control.dalpha; % differential gain
control.error = 0;
control.error_I = 0;


% create atmosphere object
atm = marsatmosphere();
% create aerocoef object
aero_coef = aeroProperties('irve');


