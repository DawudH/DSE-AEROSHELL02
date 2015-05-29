%% Definition of variables
% Points with respect to axisymmetric surface center
cg=[0.3,0,1.5]; % cg location
flap_centroid=[5, 0, 6]; % flap geometrical centroid location

% Control system limitations
phi_diff_max=deg2rad(20); % C.G. rotation system maximum rotation angle
dalpha_dt=0.2; % deg/s

% Relations between trim angle and C.G. offset, gebeund uit modnewtonian
alpha_max=25; % in degrees
cg_offset_max=1.284; % in meters, corresponds to alpha_max
CMY_A_max=-153.8; % at alpha_max

cgoffset_per_deg=cg_offset_max/alpha_max;
CMY_A_per_deg=CMY_A_max/alpha_max;

%% C.G. control section
phi_1=atan(cg(1)/cg(3)); % Default CG rotation offset angle
z_cg_max=cg(3)*tan(phi_1+phi_diff_max);
cg_alpha_trim_max=z_cg_max/cgoffset_per_deg

%% Flap control section
alpha_trim_flaps=alpha_max-cg_alpha_trim_max;
CMY_A_trim_flaps=CMY_A_per_deg*alpha_trim_flaps;
alpha_flaps_range=linspace(alpha_trim_flaps,alpha_max,10);
CMY_A_flaps_range=CMY_A_per_deg*alpha_flaps_range;
% Flap angles
theta=60; % surface inclination in degrees, w.r.t. body-symmetric axis (body frame)
alpha=cg_alpha_trim_max; % angle of attack

q=2000; % Dynamic pressure
gamma=1.4;
M=30;
num_flaps=2;

Cp_max = 2/(gamma*M^2)*((((gamma+1)^2*M^2)/(4*gamma*M^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M^2)/(gamma+1))-1);
flap_normal=[-cos(deg2rad(theta+alpha)),0,-sin(deg2rad(theta+alpha))]; % in body frame
CR=-Cp_max*(sin(deg2rad(90+theta+alpha))^2).*flap_normal;

CM_flap=cross(flap_centroid-cg,CR) % CM delivered per flap