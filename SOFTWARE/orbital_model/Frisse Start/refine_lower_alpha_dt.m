clc
clear
close all

poolobj = gcp('nocreate');
format long

%%Input constants & variables
variables

files = {'orbit_alpha_dt_isotensoid_l.txt', 'orbit_alpha_dt_apollo_l.txt', 'orbit_alpha_dt_torus_l.txt', 'orbit_alpha_dt_pastille_l.txt', 'orbit_alpha_dt_irve_l.txt'};
names = {'isotensoid', 'apollo', 'torus', 'pastille', 'irve'};

for bla = 1:length(names)
% booleans (no control, stop at beginning of elliptic orbit)
use_control = false;
multiple_orbits = false;

alpha = -25:5:25;
file_name = files{bla};

clear aero_coef
aero_coef = aeroProperties(names{bla});

%Initial Position
rx = -4100000;
ry = 10*R_m;

accuracy = 50;
init_step = 5000;
filestr = cell(length(alpha),1);
parfor i = 1:length(alpha)
    
    rx_in = rx;
    R = [rx_in,ry,0];
    go = true;
    refinestep = init_step;
    
    % init prev_c
    prev_c.t_end = false;
    prev_c.orbit = false;
    prev_c.in_atmos = false;
    prev_c.crash = false;
    prev_c.flyby = false;
    
    while go
        % run :)
        [out] = full_orbit_find(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, alpha(i)*pi/180);
        
        % save output
        filestr{i} =  [filestr{i}  sprintf(' %f %f %f %f %f %d %d %d %f \n',rx_in, ry, v, dt_atmos, alpha(i), out.c.crash, out.c.flyby, out.c.orbit, out.maxaccel)];
            
        if refinestep <= accuracy
            go = false;
            break;
        end
        if out.c.crash
            % to close
            rx_in = rx_in - refinestep;
        elseif out.c.flyby
            % to far
            rx_in = rx_in + refinestep;
        elseif out.c.orbit
            % go further to find minimum load
            rx_in = rx_in - refinestep;
        end
        
        % generate new R
        R = [rx_in,ry,0];
        
        if (prev_c.crash && out.c.flyby) || (prev_c.orbit && out.c.flyby)
            refinestep = refinestep / 2;
        end
        prev_c = out.c;
        
    end
end



%% write to file
fid = fopen(file_name,'a');
for k = 1:length(alpha)
    fprintf(fid,'%s',filestr{k});
end
fclose(fid);

end
% close parallel processing
delete(poolobj);
