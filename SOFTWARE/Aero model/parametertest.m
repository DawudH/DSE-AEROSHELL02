clear; clc; close all;


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
%flat:
% x = [0.0    0.01   -0.1556    0.5029    2.6577   -0.5685    0.3564    2.1639   -0.7676    2.4901];

%rounder:
x = [0   3   -9.3661   -8.5449   31.2249  -23.2111   -6.7477    9.3702   25.7382    0.0021];
radiusglobal = 1;
LoverDglobal = -0.3;
qglobal = 30;

heightfactor = x(2);
height = radiusglobal * heightfactor;

skewness = 0; %x(1);

poly = [x(3:end),0,0];

gamma = 1.29; %gamma
rho = 1e-5; %Density
T = 150; %Temperatuure
beta = 0; %Sideslip
q = 30;

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
beta0data.CRAX =  mod.CRA_body_array(1,:);
beta0data.CRAY =  mod.CRA_body_array(2,:);
beta0data.CRAZ =  mod.CRA_body_array(3,:);
beta0data.CRAD =  mod.CRA_aero_array(1,:);
beta0data.CRAS =  mod.CRA_aero_array(2,:);
beta0data.CRAL =  mod.CRA_aero_array(3,:);
beta0data.CMAX =  mod.CMA_body_array(1,:);
beta0data.CMAY =  mod.CMA_body_array(2,:);
beta0data.CMAZ =  mod.CMA_body_array(3,:);
beta0data.CRX =  mod.CR_body_array(1,:);
beta0data.CRY =  mod.CR_body_array(2,:);
beta0data.CRZ =  mod.CR_body_array(3,:);
beta0data.CRD =  mod.CR_aero_array(1,:);
beta0data.CRS =  mod.CR_aero_array(2,:);
beta0data.CRL =  mod.CR_aero_array(3,:);
beta0data.CMX =  mod.CM_body_array(1,:);
beta0data.CMY =  mod.CM_body_array(2,:);
beta0data.CMZ =  mod.CM_body_array(3,:);

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CM_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CM_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CM_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_M [-]$', 'interpreter', 'latex');
legend('C_{M_x} [-]', 'C_{M_y} [-]', 'C_{M_z} [-]', 'interpreter', 'latex');
title('Moment, Beta=0deg');

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CR_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CR_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CR_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{R_body} [-]$', 'interpreter', 'latex');
legend('C_{X} [-]', 'C_{Y} [-]', 'C_{Z} [-]', 'interpreter', 'latex');
title('Force, body frame, beta=0deg');

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{R_aerodynamic} [-]$', 'interpreter', 'latex');
legend('C_{D} [-]', 'C_{S} [-]', 'C_{L} [-]', 'interpreter', 'latex');
title('Force, aero frame, beta=0deg');


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
beta30data.CRAX =  mod.CRA_body_array(1,:);
beta30data.CRAY =  mod.CRA_body_array(2,:);
beta30data.CRAZ =  mod.CRA_body_array(3,:);
beta30data.CRAD =  mod.CRA_aero_array(1,:);
beta30data.CRAS =  mod.CRA_aero_array(2,:);
beta30data.CRAL =  mod.CRA_aero_array(3,:);
beta30data.CMAX =  mod.CMA_body_array(1,:);
beta30data.CMAY =  mod.CMA_body_array(2,:);
beta30data.CMAZ =  mod.CMA_body_array(3,:);
beta30data.CRX =  mod.CR_body_array(1,:);
beta30data.CRY =  mod.CR_body_array(2,:);
beta30data.CRZ =  mod.CR_body_array(3,:);
beta30data.CRD =  mod.CR_aero_array(1,:);
beta30data.CRS =  mod.CR_aero_array(2,:);
beta30data.CRL =  mod.CR_aero_array(3,:);
beta30data.CMX =  mod.CM_body_array(1,:);
beta30data.CMY =  mod.CM_body_array(2,:);
beta30data.CMZ =  mod.CM_body_array(3,:);

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CM_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CM_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CM_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_M [-]$', 'interpreter', 'latex');
legend('C_{M_x} [-]', 'C_{M_y} [-]', 'C_{M_z} [-]', 'interpreter', 'latex');
title('Moment, Beta=30deg');

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CR_body_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CR_body_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CR_body_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{R_body} [-]$', 'interpreter', 'latex');
legend('C_{X} [-]', 'C_{Y} [-]', 'C_{Z} [-]', 'interpreter', 'latex');
title('Force, body frame, beta=30deg');

figure;
hold on;
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(1,:));
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(2,:));
plot(rad2deg(mod.alpha_array), mod.CR_aero_array(3,:));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{R_aerodynamic} [-]$', 'interpreter', 'latex');
legend('C_{D} [-]', 'C_{S} [-]', 'C_{L} [-]', 'interpreter', 'latex');
title('Force, aero frame, beta=30deg');

save('outputfiles/stability_round.mat', 'x', 'mod', 'beta0data', 'beta30data');

mod.geom.plotGeometry(true, false);
