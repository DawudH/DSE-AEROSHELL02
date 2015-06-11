clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = false;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

control.alpha_init = 15/180*pi;
gamma = 21.845;
i=1;

for k = 0:2:180
    %%function
    [out] = full_orbit(k,R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
    if out.c.orbit
        results.mu(i) = k;
        results.h(i) = min(sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m);
        results.q(i) = max(out.q);
        results.M(i) = max(out.M);
        i = i + 1;
    end
end
save('effect_mu.mat','results')