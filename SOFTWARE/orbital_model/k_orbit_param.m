function [kepler_param] = k_orbit_param(R1,R2,V,dt,G,M)

%% Determine kepler orbit parameters
mu = G*M;
% get magintudes
r = norm(R1);
v = norm(V);
r2 = norm(R2);

% calculate the semi major axis with visviva
a = r*mu / (2*mu - r * v^2);

% determine e with areal velocity
dtheta = atan2(R2(2),R2(1)) - atan2(R1(2),R1(1))
dA = 1/2 * r * r2 * sin(dtheta)
e = sqrt( 1 - 4 * (dA / dt)^2 / (a * mu) );

% calculate theta
theta = acos( (a*(1-e^2)-r) / (e*r) );

if a>0
    % calculate the orbital period
    T = 2*pi * sqrt(a^3 / mu);

    % semi-minor axis: b
    b = a*sqrt(1-e^2);
else
    % calculate the orbital period
    T = Inf;

    % semi-minor axis: b
    b = a*sqrt(e^2-1);
end
% perifocal distance
rp  = a*(1-e);

% apofocal distance
ra = 2*a - rp;

% 
theta0 = atan2(R1(2),R1(1));
thetap = theta0 - theta;
    
% output
kepler_param.e = e; %
kepler_param.a = a; % 
kepler_param.b = b; %
kepler_param.theta = theta; %
kepler_param.thetap = thetap; %
kepler_param.theta0 = theta0; %
kepler_param.rp = rp; %
kepler_param.ra = ra; %
kepler_param.T = T; %

kepler_param
end