clc
clear all
close all

% load constants
constants

ry = 10* R_m; %[m]
v = 5000; %[m/s]
dt_init = 1;
dt_atmos = 1;
dt_kep_init = 1e-8;
CD = 1.25;
rx = -4682000.0;
CL = 0.56;
tend = 3600 * 3;

[out] = orbit_full(rx,ry,CD,v,dt_init,dt_atmos,dt_kep_init,CL,tend);