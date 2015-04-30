function [out] = orbit_selection(rx,ry,CD,v,dt,CL,R_m,Omega_m,S,m,G,M_mars,h_atm,crash_margin,g_earth)


% define output 
out.inorbit = false;
out.inatmos = false;
out.crash = false;

% boundary conditions
out.R(1,:) = [rx, ry, 0];
out.V(1,:) = [0, -v, 0];
out.a(1,:) = [0, 0, 0];
out.ag(1,:) = [0, 0, 0];
out.ad(1,:) = [0, 0, 0];
out.al(1,:) = [0, 0, 0];

% create atmosphere object
atm = marsatmosphere();

% escape velocity at the border of the atmosphere
V_esc = sqrt(G*M_mars * 2 / (h_atm + R_m)); 

% start while loop until crash, orbit or passby..
i = 1;
while true
    
    % get orbital parameters at next node
    orbit_new = orbit(out.R(i,:),out.V(i,:),out.a(i,:),CD,CL,dt,atm,R_m,Omega_m,S,m);
    out.R(i+1,:) = orbit_new.R;
    out.V(i+1,:) = orbit_new.V;
    out.a(i+1,:) = orbit_new.a;
    out.ad(i+1,:) = orbit_new.ad;
    out.al(i+1,:) = orbit_new.al;
    out.ag(i+1,:) = orbit_new.ag;
        
    % if the s/c has ever been in the atmosphere in_atmus => true
    if (norm(out.R(i+1,:)) < (R_m + h_atm)) 
        out.inatmos = true;
    end
        
    % if s/c has been in atmos and is leaving atmos again (at next node), stop simulation
    % and check if the exit velocity is smaller then the escape velocity
    % (then it is in orbit!)
    if out.inatmos && ((norm(out.R(i+1,:)) > (R_m + h_atm))) 
        
        if (norm(out.V(i+1,:)) < V_esc)
            out.inorbit = true; 
        end
        break;
        
    end
    
    % check if crashed.. % crashed = R_m + crash margin (due to to big timesteps) 
     if (out.inatmos) && (norm(out.R(i+1,:)) < (R_m + crash_margin) ) 
        out.crash = true;
        break;
    end

    % if passed the planet without getting into the atmosphere..
    if (out.inatmos == false) && (norm(out.R(i+1,2)) < -(R_m + h_atm))
        % escape velocity at the inital point
        V_esc_0 = sqrt(G*M_mars * 2 / norm(out.R(1,:))); 
        % check if initial velocity is sufficient to get into orbit
        if v < V_esc_0
            out.inorbit = true;
        end
        break;
    end
    
    
    % counter :)
    i = i+1;
end

% calculate the magntude of the acceleration experienced by the humans
out.a_human_mag = sqrt((out.ad(:,1)+out.al(:,1)).^2 + (out.ad(:,2)+out.al(:,2)).^2 + (out.ad(:,3)+out.al(:,3)).^2);
out.maxaccel = max(out.a_human_mag)/g_earth;
disp(['rx = ' num2str(rx) ' [m], CD = ' num2str(CD) ' [-], CL = ' num2str(CL) ' [-], in atmosphere: ' num2str(out.inatmos) ', crashed: ' num2str(out.crash) ', in orbit: ' num2str(out.inorbit) ', acceleration: ' num2str(out.maxaccel) ])
end
