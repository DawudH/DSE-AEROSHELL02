clc
clear all
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = true;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

%%function
%[out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
load('../../Aero model/orbits/orbit_iteration_0_1.mat')
mu = G*M_mars;
r = out.okep.a*(1-out.okep.e^2)/(1-out.okep.e)
v = sqrt(mu*(2/r-1/out.okep.a))
a_n = (R_m+h_atm+r)/2;
v_n = sqrt(mu*(2/r-1/a_n))
%v_c = sqrt(mu/r)
%dv_c = v_c-v
%v_n_p = sqrt(mu*(2/(R_m+h_atm)-1/a_n))
%v_c_p = sqrt(mu/(R_m+h_atm))
dv1 = v_n - v
%dv2 = v_c_p - v_n_p