clc
clear all
close all

% load constants
constants

ry = 10* R_m; %[m]
v = 6000; %[m/s]
dt_init = 1;
dt_atmos = 0.1;
dt_kep_init = 1e-8;
CD = 1.25;
rx = -4383333.333333;
CL = -0.062500 ;
tend = 3600 * 24* 10;

[out] = orbit_full(rx,ry,CD,v,dt_init,dt_atmos,dt_kep_init,CL,tend);