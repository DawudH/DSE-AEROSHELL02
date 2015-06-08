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


alpha = -10:2.5:20; 
gamma_range = 21.8:0.01:24;
k = 1;
for j = 1:length(alpha)
    for i = 1:length(gamma_range)
        
        control.alpha_init = alpha(j)*pi/180; % rad
        gamma = gamma_range(i);
        
        %%function
        [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);

        if out.c.orbit
            results.ALPHA(k) = control.alpha_init;
            results.GAMMA(k) = gamma;
            results.A(k) = norm(max(out.A_aero));
            k = k+1;
        elseif  out.c.crash
            break;
        end
    end
end

save('results.mat','results')
plot(results.GAMMA,results.ALPHA*180/pi,'x')


    
