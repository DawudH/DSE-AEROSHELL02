clc
clear
close all

%%Input constants & variables
variables
step = -100;
rx_b = -4144000;
rx_e = -4200000;
rx = (rx_b:step:rx_e)'; %[m]
n = length(rx);
ry = 10*R_m*ones(n,1); %[m]
R = [rx,ry,zeros(n,1)];
out_refine.rx = rx(1);
out_orbit.c.orbit = false;

% booleans
use_control = true;
multiple_orbits = false;

while (out_orbit.c.orbit == false)
    for i=1:length(rx)
    %%function
        [out_orbit] = full_orbit(R(i,:), V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef,use_control,multiple_orbits);
        i
        out_orbit.c
        [out_refine] = refine_check(out_orbit.c,rx(i),out_refine);
        if out_refine.refine
            step = (rx(i)-out_refine.rx)/10;
            rx = (out_refine.rx:step:rx(i))';
            R = [rx,ry,zeros(n,1)];
            break;
        end
        if out_orbit.c.orbit
            rx(i)
        end
        if out_refine.break
            break;
        end
    end
    if out_refine.refine == false
        disp('no orbit found')
        break;
    end
end