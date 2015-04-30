function [t, conditions] = full_orbit()
constants
conditions.atmos = false; conditions.crash = false; conditions.orbit = true; 
rx = -4159317; %[m]
ry = 10*R_m; %[m]
v = 7000; %[m/s]
R = [rx,ry,0]; %[m]
V = [0,-v,0]; %[m/s]
a = [0,0,0]; %[m/s^2]
CD = 1.2; %[-]
CLCD = 0.3; %[-]
CL = CLCD*CD; %[-]
dt = 1; %[s]
new = orbit(R,V,a,CD,CL,dt);
V_esc = sqrt(G*M_mars * 2 / (h_atm + R_m)); % escape velocity at the border of the atmosphere
V_esc0 = sqrt(G*M_mars * 2 / norm(R));
i = 1;
while true
    %calculate new orbital parameters
    new = orbit(new.R,new.V,new.a,CD,CL,dt);
    
    %check for being in the atmosphere
    if norm(new.R)<(R_m+h_atm)
        conditions.atmos = true;
        %check for crash
        if norm(new.R)<R_m
            conditions.crash = true;
            t = i*dt;
            break;
        end
    end
    
    %check if s/c has been in atmosphere and comes out faster than the escape velocity
    if (norm(new.R)>(R_m+h_atm)) && (atmos == true) && (norm(V)>V_esc)
        conditions.orbit = false;
        t = i*dt;
        break;
    end
    
    %check if s/c didn't go past mars without touching the atmosphere
    if (new.R(2)<-R_m) && (atmos == false) && (v>V_esc0)
        conditions.orbit = false;
        t = i*dt;
        break;
    end
    
    
    i = i+1