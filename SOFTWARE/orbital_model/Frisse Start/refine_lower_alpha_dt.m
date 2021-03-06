clc
clear
close all

poolobj = gcp('nocreate');
format long

%%Input constants & variables
variables

%files = {'orbit_alpha_dt_isotensoid_l.txt', 'orbit_alpha_dt_apollo_l.txt', 'orbit_alpha_dt_torus_l.txt', 'orbit_alpha_dt_pastille_l.txt', 'orbit_alpha_dt_irve_l.txt'};
%names = {'isotensoid', 'apollo', 'torus', 'pastille', 'irve'};

files = {'orbit_alpha_dt_irve_alpha_20.txt'};
names = {'irve'};

for bla = 1:length(names)
% booleans (use control, stop at beginning of elliptic orbit)
use_control = true;
multiple_orbits = false;

alpha = 20; % deg
dalphadt = -1.5:0.25:-0.5; %deg
file_name = files{bla};

clear aero_coef
aero_coef = aeroProperties(names{bla});

%Initial Position
rx = -4150000;
ry = 10*R_m;

accuracy = 10;
init_step = 500;

filestr = cell(length(dalphadt),1);
parfor i = 1:length(dalphadt)
    
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
        [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, alpha*pi/180,dalphadt(i)*pi/180);

        % save output
        filestr{i} =  [filestr{i}  sprintf(' %f %f %f %f %f %d %d %d %f %f %f %f %f \n',rx_in, ry, v, dt_atmos, dalphadt(i), out.c.crash, out.c.flyby, out.c.orbit, out.maxaccel, out.alpha(1), out.Kp, out.Ki, out.Kd)];
  
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
        
        if (prev_c.crash && out.c.flyby) || (prev_c.orbit && out.c.flyby) || (prev_c.flyby && out.c.orbit) || (prev_c.flyby && out.c.crash)
            refinestep = refinestep / 2;
        end
        prev_c = out.c;
        
    end
end



%% write to file
fid = fopen(file_name,'a');
for k = 1:length(dalphadt)
    fprintf(fid,'%s',filestr{k});
end
fclose(fid);

end
% close parallel processing
delete(poolobj);
