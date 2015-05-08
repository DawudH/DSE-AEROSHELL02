function [out] = orbit_full(rx,ry,CD,v,dt_init,dt_atmos,dt_kep_init,CL,tend);
% load constants
constants

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
% waitbar
h = waitbar(0,'Initializing waitbar...');

% kepler orbit 
to_kepler = false;

% define time step
dt = dt_init;

while true
    
    if to_kepler
        
        [orbit_new,t_kep] = orbit_kepler(kepler_param,orbit_new);
        to_kepler = false;
        out.inorbit = false;
        
    else
        % get orbital parameters at next node
        orbit_new = orbit(out.R(i,:),out.V(i,:),out.a(i,:),CD,CL,dt,atm,R_m,Omega_m,S,m);
    end
    
    out.R(i+1,:) = orbit_new.R;
    out.V(i+1,:) = orbit_new.V;
    out.a(i+1,:) = orbit_new.a;
    out.ad(i+1,:) = orbit_new.ad;
    out.al(i+1,:) = orbit_new.al;
    out.ag(i+1,:) = orbit_new.ag;

    if out.inorbit && (to_kepler == false)
        
        orbit_new = orbit(out.R(i,:),out.V(i,:),out.a(i,:),CD,CL,dt_kep_init,atm,R_m,Omega_m,S,m);
        
        dtheta = atan2(orbit_new.R(2),orbit_new.R(1)) - atan2(out.R(i,2),out.R(i,1));
        kepler_param = k_orbit_param(out.R(i,:),orbit_new.R,out.V(i,:),dt_kep_init,dtheta,G,M_mars);
        % BEUN stop om te controleren of script klopt
        disp(kepler_param);
        to_kepler = true;
        out.inatmos = false;
        % break;
    end
        
        
        
    % if the s/c has ever been in the atmosphere in_atmus => true
    % change time step to time step in atmos!
    if (norm(out.R(i+1,:)) < (R_m + h_atm)) 
        out.inatmos = true;
        dt = dt_atmos;
    end
        
    % if s/c has been in atmos and is leaving atmos again (at next node)
    % check if the exit velocity is smaller then the escape velocity
    % (then it is in orbit!) otherways end the simulation
    if (out.inatmos) && ((norm(out.R(i+1,:)) > (R_m + h_atm))) 
        
        if (norm(out.V(i+1,:)) < V_esc)
            out.inorbit = true; 
        else
            break;
        end
        
        
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
        else
            break;
        end
        
    end
    
    % end of simulation time!
    if (i*dt >= tend)
    
        break;
        
    end
    
    
    % counter :)
    i = i+1;
    
    
    %progress bar
    perc = i*dt / tend * 100;
    waitbar(perc/100,h,sprintf('%5.1f %...',perc))
end

%close waitbar
close(h);

% calculate the magntude of the acceleration experienced by the humans
out.a_human_mag = sqrt((out.ad(:,1)+out.al(:,1)).^2 + (out.ad(:,2)+out.al(:,2)).^2 + (out.ad(:,3)+out.al(:,3)).^2);
out.maxaccel = max(out.a_human_mag)/g_earth;
disp(['rx = ' num2str(rx) ' [m], CD = ' num2str(CD) ' [-], CL = ' num2str(CL) ' [-], in atmosphere: ' num2str(out.inatmos) ', crashed: ' num2str(out.crash) ', in orbit: ' num2str(out.inorbit) ', acceleration: ' num2str(out.maxaccel) ])

% plot orbit
% circle plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
figure('name','Orbit')
grid on
axis equal
hold on
plot(out.R(:,1),out.R(:,2))
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g')
% plot kepler orbit
if exist('kepler_param','var')
   
    rk = kepler_param.a * (1- kepler_param.e^2) ./ (1 + kepler_param.e .* cos(theta_plot));
    polar(theta_plot+kepler_param.thetap,rk,'k');
    plot(kepler_param.rp*cos(kepler_param.thetap),kepler_param.rp*sin(kepler_param.thetap),'*')
    plot(-kepler_param.ra*cos(kepler_param.thetap),-kepler_param.ra*sin(kepler_param.thetap),'d')
    plot((R_m + h_atm)*cos(kepler_param.thetap+kepler_param.theta),(R_m + h_atm)*sin(kepler_param.thetap+kepler_param.theta),'p')
    plot((R_m + h_atm)*cos(kepler_param.thetap-kepler_param.theta),(R_m + h_atm)*sin(kepler_param.thetap-kepler_param.theta),'p')
end


t = 0:dt:(length(out.R)*dt - dt);
figure('name','parameters over time')
subplot(3,1,1)
Rm = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2);
plot(t,Rm)
grid on
subplot(3,1,2)
Vm = sqrt(out.V(:,1).^2 + out.V(:,3).^2 + out.V(:,2).^2);
plot(t,Vm)
grid on
subplot(3,1,3)
am = sqrt(out.a(:,1).^2 + out.a(:,3).^2 + out.a(:,2).^2);
plot(t,am)
grid on
end
