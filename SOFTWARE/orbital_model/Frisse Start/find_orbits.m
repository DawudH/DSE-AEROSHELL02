clc
clear
close all

%%Input constants & variables
variables
step = -100;
rx_b = -4720500;
rx_e = -4725000;
rx = (rx_b:step:rx_e)'; %[m]
n = length(rx);
ry = SOI*ones(n,1); %[m]
R = [rx,ry,zeros(n,1)];
out_refine.rx = rx(1);
out_orbit.c.orbit = false;

while (out_orbit.c.orbit == false) && (abs(step)>1/10000)
    for i=1:length(rx)
    %%function
        [out_orbit] = full_orbit(R(i,:), V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth);
        i
        out_orbit.c
        [out_refine] = refine_check(out_orbit.c,rx(i),out_refine);
        if out_refine.refine
            step = (rx(i)-out_refine.rx)/10;
            rx = (out_refine.rx:step:rx(i))';
            n = length(rx);
            ry = SOI*ones(n,1); %[m]
            R = [rx,ry,zeros(n,1)];
            break;
        end
        if out_orbit.c.orbit
            rx(i)
        end
    end
end
for i=1:length(rx)
%%function
    [out_orbit] = full_orbit(R(i,:), V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth);
    out_orbit.c
    [out_refine] = refine_check(out_orbit.c,rx(i),out_refine);
    if out_refine.refine
        break;
    end
    if out_orbit.c.orbit
        out_orbit.c
    end
end