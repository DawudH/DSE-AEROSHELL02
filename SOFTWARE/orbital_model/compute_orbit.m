clc
clear all
close all

% load constants
constants

ry = SOI; %[m]
v = 7000; %[m/s]
dt_init = 1;
dt_atmos = 0.5;
dt_kep_init = 1e-8;
CD = 1.25;
rx =  -4389523.809524;
CL = -0.187500 ;
tend = 3600 * 24* 5;

[out] = orbit_full(rx,ry,CD,v,dt_init,dt_atmos,dt_kep_init,CL,tend);