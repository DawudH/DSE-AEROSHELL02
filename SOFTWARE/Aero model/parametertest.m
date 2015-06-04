clear; clc;


%% Get the minimum and maximum dynamic pressure
d12_one_time.qmin = inf;
d12_one_time.qmax = -inf;
d12_one_time.Mmin = inf;
d12_one_time.Mmax = -inf;
d12_one_time.alphamin = inf;
d12_one_time.alphamax = -inf;

d6_one_time.qmin = inf;
d6_one_time.qmax = -inf;
d6_one_time.Mmin = inf;
d6_one_time.Mmax = -inf;
d6_one_time.alphamin = inf;
d6_one_time.alphamax = -inf;

d12_just_orbit.qmin = inf;
d12_just_orbit.qmax = -inf;
d12_just_orbit.Mmin = inf;
d12_just_orbit.Mmax = -inf;
d12_just_orbit.alphamin = inf;
d12_just_orbit.alphamax = -inf;

d6_just_orbit.qmin = inf;
d6_just_orbit.qmax = -inf;
d6_just_orbit.Mmin = inf;
d6_just_orbit.Mmax = -inf;
d6_just_orbit.alphamin = inf;
d6_just_orbit.alphamax = -inf;

load('orbits\out_d12_one_time.mat');
d12_one_time.qmin = min(d12_one_time.qmin, min(out.q));
d12_one_time.qmax = max(d12_one_time.qmax, max(out.q));
d12_one_time.Mmin = min(d12_one_time.Mmin, min(out.M));
d12_one_time.Mmax = max(d12_one_time.Mmax, max(out.M));
d12_one_time.alphamin = rad2deg(min(d12_one_time.alphamin, min(out.alpha)));
d12_one_time.alphamax = rad2deg(max(d12_one_time.alphamax, max(out.alpha)));

load('orbits\out_d6_one_time.mat');
d6_one_time.qmin = min(d6_one_time.qmin, min(out.q));
d6_one_time.qmax = max(d6_one_time.qmax, max(out.q));
d6_one_time.Mmin = min(d6_one_time.Mmin, min(out.M));
d6_one_time.Mmax = max(d6_one_time.Mmax, max(out.M));
d6_one_time.alphamin = rad2deg(min(d6_one_time.alphamin, min(out.alpha)));
d6_one_time.alphamax = rad2deg(max(d6_one_time.alphamax, max(out.alpha)));

load('orbits\out_d12_just_orbit.mat');
d12_just_orbit.qmin = min(d12_just_orbit.qmin, min(out.q));
d12_just_orbit.qmax = max(d12_just_orbit.qmax, max(out.q));
d12_just_orbit.Mmin = min(d12_just_orbit.Mmin, min(out.M));
d12_just_orbit.Mmax = max(d12_just_orbit.Mmax, max(out.M));
d12_just_orbit.alphamin = rad2deg(min(d12_just_orbit.alphamin, min(out.alpha)));
d12_just_orbit.alphamax = rad2deg(max(d12_just_orbit.alphamax, max(out.alpha)));

load('orbits\out_d6_just_orbit.mat');
d6_just_orbit.qmin = min(d6_just_orbit.qmin, min(out.q));
d6_just_orbit.qmax = max(d6_just_orbit.qmax, max(out.q));
d6_just_orbit.Mmin = min(d6_just_orbit.Mmin, min(out.M));
d6_just_orbit.Mmax = max(d6_just_orbit.Mmax, max(out.M));
d6_just_orbit.alphamin = rad2deg(min(d6_just_orbit.alphamin, min(out.alpha)));
d6_just_orbit.alphamax = rad2deg(max(d6_just_orbit.alphamax, max(out.alpha)));
d12_one_time
d6_one_time
d12_just_orbit
d6_just_orbit

clear;
x = [4.8627   10.8197   22.6692   -9.3661   -8.5449   31.2249  -23.2111   -6.7477    9.3702   25.7382    0.0021];
radiusglobal = 1;
LoverDglobal = -0.3;
qglobal = 30;

heightfactor = 3;%x(2);
height = radiusglobal * heightfactor;

skewness = 0; %x(1);

poly = [x(3:end),0,0];

gamma = 1.29; %gamma
rho = 1e-5; %Density
T = 150; %Temperatuure
q = 40; % Maat voor aantal elementen

beta = 0; %Sideslip

phi = 0;
V = 7500;
a = 187.5;
radius = 6;
center = [0 0 0];
[TriGeom, A] = ParaGeom(q, skewness, radius, height, poly);
geom = aeroGeometry(TriGeom, A);
mod = modnewtonian(geom, gamma, a, center, rho, T);

%% Angle of attack sweep, beta=0 deg
% Angle of attack values
alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 40; %degrees
beta = 0;
phi = 0;

mod = mod.alphasweep(V, beta, phi, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

beta0data.alpha = mod.alpha_array;
beta0data.CMAX =  mod.CMA_body_array(1,:);
beta0data.CMAY =  mod.CMA_body_array(2,:);
beta0data.CMAZ =  mod.CMA_body_array(3,:);

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_m [-]$', 'interpreter', 'latex');
legend('C_{m_x} [-]', 'C_{m_y} [-]', 'C_{m_z} [-]', 'interpreter', 'latex');
title('beta=0deg');



%% Angle of attack sweep, beta=30 deg
mod = modnewtonian(geom, gamma, a, center, rho, T);
% Angle of attack values
alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 40; %degrees
beta = deg2rad(30);
phi = 0;

mod = mod.alphasweep(V, beta, phi, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

beta30data.alpha = mod.alpha_array;
beta30data.CMAX =  mod.CMA_body_array(1,:);
beta30data.CMAY =  mod.CMA_body_array(2,:);
beta30data.CMAZ =  mod.CMA_body_array(3,:);


figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CMA_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_m [-]$', 'interpreter', 'latex');
legend('C_{m_x} [-]', 'C_{m_y} [-]', 'C_{m_z} [-]', 'interpreter', 'latex');
title('beta=30deg');

save('outputfiles/stability.mat', 'x', 'mod', 'beta0data', 'beta30data');
