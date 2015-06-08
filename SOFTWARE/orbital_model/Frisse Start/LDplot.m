clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = true;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

CLCD = -1:0.01:0;
aero_coef.cda = 138.5;
figure('name','LD')
hold on
for k = 1:length(CLCD)
    aero_coef.cla = aero_coef.cda*CLCD(k);
    [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
    results.CLCD(k) = CLCD(k);
    results.A(k) = norm(max(out.A_aero));
    results.orbit(k) = out.c.orbit;
    if results.orbit(k)
        plot(results.CLCD(k),results.A(k)/g_earth,'bx')
    else
        plot(results.CLCD(k),results.A(k)/g_earth,'rx')
    end
end
save('resultsLD60.mat','results')
