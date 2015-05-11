clc
clear all
close all

% load constants
constants

ry = 10* R_m; %[m]
v = 6000; %[m/s]
dt_init = 1;
dt_atmos = 0.5;
dt_kep_init = 1e-8;
rx = -4386666;
tend = 3600 * 24* 1;

control.CL_range = [-0.35 0.35];
control.CLCD = 0.25;
control.a = 2.9*g_earth;
control.CLa = 0.02;
control.dalpha = 0.2;
control.CL_init = -0.18;



[out] = orbit_full(rx,ry,v,dt_init,dt_atmos,dt_kep_init,tend,control);