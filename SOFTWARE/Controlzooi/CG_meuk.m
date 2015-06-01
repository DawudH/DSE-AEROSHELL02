%% Definition of variables
% Points with respect to axisymmetric surface center
cg=[0.3,0,3.]; % cg location
flap_centroid=[5, 0, 6]; % flap geometrical centroid location

% Control system limitations
phi_diff_max=deg2rad(10); % C.G. rotation system maximum rotation angle
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
cg(1)=z_cg_max;
cg_alpha_trim_max=z_cg_max/cgoffset_per_deg

%% Flap control section
alpha_trim_flaps=alpha_max-cg_alpha_trim_max; % alpha_trim handled by flaps
CMY_A_trim_flaps=CMY_A_per_deg*alpha_trim_flaps; % CMY_A handled by flaps to trim to alpha_max
%alpha_flaps_range=linspace(cg_alpha_trim_max,alpha_max,3);
%CMY_A_flaps_range=CMY_A_per_deg*alpha_flaps_range;

flap_offset=flap_centroid-cg;

% Flap anglescg_alpha_trim_max
%theta=60; % surface inclination in degrees, w.r.t. body-symmetric axis (body frame)
theta_range=linspace(10,60,4);

q=2000; % Dynamic pressure
gamma=1.4;
M=30;

Cp_max = 2/(gamma*M^2)*((((gamma+1)^2*M^2)/(4*gamma*M^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M^2)/(gamma+1))-1);
flap_normal=[-cos(deg2rad(theta_range+alpha_max))',zeros(size(theta_range))',-sin(deg2rad(theta_range+alpha_max))'];
CR=flap_normal'*diag(-Cp_max*(sin(deg2rad(90+alpha_max+theta_range)).^2));

CM_flap=cross(repmat(flap_offset',1,length(theta_range)),CR); % CM delivered per flap

plot(theta_range,CM_flap(2,:)');