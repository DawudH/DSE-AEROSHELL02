clc
clear
close all

%%Input constants & variables
variables
%escape velocity at boundary of atmosphere
V_esc = sqrt(2*G*M_mars/(R_m+h_atm));

%%function
%[out] = full_orbit(R, V, V_esc, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control);
                    
%%processing (plot/write to file)


hyperbolic_kepler(R, V,A,G,M_mars,R_m,h_atm,dt_kep_init)