clc
clear all
close all

%%Input constants & variables
variables
%escape velocity at boundary of atmosphere
V_esc = sqrt(2*G*M/(R_m+h_atm));

%%function
[out] = full_orbit(R,V,A,G,M_mars,R_m,h_atm,dt_kep_init,dt_atmos);

%%processing (plot/write to file)
