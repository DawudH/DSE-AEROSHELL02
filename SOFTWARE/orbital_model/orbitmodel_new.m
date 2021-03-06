function [out, R, V, a] = orbitmodel_new(rx,ry,R_m,m,CD,S,v,dt,h_atmos,M_mars,G,CLCD)


omega_m = 2*pi / (24.6229622 * 3600 ); %[rad/sec]
Omega_m = [0,0,1]*omega_m;

% define output 
out.inorbit = false;
out.inatmos = false;
out.maxaccel = 0;
out.crash = false;

R(1,:) = [rx, ry, 0];
V(1,:) = [0, -v, 0];
Vatm(1,:) = cross(Omega_m,R(1,:));
[ g, p, T, rho, asound ] = mars_atmosphere(norm(R(1,:)) - R_m);
ag(1,:) = -g * R(1,:)/norm(R(1,:));
vel_unit(1,:) = (V(1,:) - Vatm(1,:)) / norm( V(1,:) - Vatm(1,:));
ad(1,:) = - vel_unit(1,:) * CD * 0.5* rho * norm(V(1,:) - Vatm(1,:))^2 * S / m;
al(1,:) = - cross(vel_unit(1,:),[0 0 1]) * CLCD*CD * 0.5* rho * norm(V(1,:) - Vatm(1,:))^2 * S / m;
a(1,:) = ag(1,:) + ad(1,:) + al(1,:);
i = 1;
V_esc = sqrt(G*M_mars * 2 / (h_atmos + R_m)); % escape velocity at the border of the atmosphere

while true
    
    [ g, p, T, rho, asound ] = mars_atmosphere(norm(R(i,:)) - R_m);
    
    ag(i+1,:) = -g * R(i,:)/norm(R(i,:));
    ad(i+1,:) = - vel_unit(i,:) * CD * 0.5* rho * norm(V(i,:) - Vatm(i,:))^2 * S / m;
    al(i+1,:) = - cross(vel_unit(i,:),[0 0 1]) * CLCD*CD * 0.5* rho * norm(V(i,:) - Vatm(i,:))^2 * S / m;
    a(i+1,:) = ag(i+1,:) + ad(i+1,:) + al(i+1,:);
    R(i+1,:) = R(i,:) + V(i,:)*dt+a(i,:)*dt^2;
    Vatm(i+1,:) = cross(Omega_m,R(i,:));
    V(i+1,:) = V(i,:) + a(i,:)*dt;
    vel_unit(i+1,:) = (V(i+1,:) - Vatm(i+1,:)) / norm( V(i+1,:) - Vatm(i+1,:));
        
    % if the s/c has ever been in the atmosphere in_atmus => true
    if (norm(R(i,:)) < (R_m + h_atmos)) 
        out.inatmos = true;
        
    end
        
    % if s/c has been in atmos and is leaving atmos again, stop simulation
    % and check velocity and max deceleration...
    if out.inatmos && ((norm(R(i,:)) > (R_m + h_atmos))) 

        
        %disp(['Vesc = ' num2str(V_esc) ' [m/s], Vend = ' num2str(norm(V(end,:))) ' [m/s]'])
        if (norm(V(end,:)) < V_esc)
            out.inorbit = true; 
        end
        break;
        
    end
    
    % check if crashed.. % crashed = R_m + 5 km..
    if out.inatmos && (norm(R(i,:)) < R_m + 5000) 
        out.crash = true;
        break;
    end

    % if passed the planet without getting into the atmosphere..
    if (out.inatmos == false) && (R(i,2) < -(R_m + h_atmos))
        V_esc_0 = sqrt(G*M_mars * 2 / norm(R(1,:))); % escape velocity at the inital point
        % check if initial velocity is to fast to get into orbit.. ?
        if v < V_esc_0
            out.inorbit = true;
        end
        break;
    end
    
    
     i = i+1;
end
admag = sqrt((ad(:,1)+al(:,1)).^2 + (ad(:,2)+al(:,2)).^2 + (ad(:,3)+al(:,3)).^2);
out.maxaccel = max(admag)/9.81;
disp(['rx = ' num2str(rx) ' [m], CD = ' num2str(CD) ' [-], CL/CD = ' num2str(CLCD) ' [-], in atmosphere: ' num2str(out.inatmos) ', crashed: ' num2str(out.crash) ', in orbit: ' num2str(out.inorbit) ', acceleration: ' num2str(out.maxaccel) ])
end
