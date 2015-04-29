function [out, R, V, a] = orbitmodel_new(rx,ry,R_m,m,CD,S,v,dt,h_atmos,M_mars,G)


% define output 
out.inorbit = false;
out.inatmos = false;
out.rightaccel = false;
out.maxaccel = 0;
out.crash = false;
out.deltav = 0;

R(1,:) = [rx, ry, 0];
V(1,:) = [0, -v, 0];
[ g, p, T, rho, asound ] = mars_atmosphere(norm(R(1,:)) - R_m);
ag(1,:) = -g * R(1,:)/norm(R(1,:));
ad(1,:) = - (V(1,:) / norm(V(1,:))) * CD * 0.5* rho * norm(V(1,:))^2 * S / m;
a(1,:) = ag(1,:) + ad(1,:);
i = 1;
V_esc = sqrt(G*M_mars * 2 / (h_atmos + R_m)); % escape velocity at the border of the atmosphere
while true
    
    [ g, p, T, rho, asound ] = mars_atmosphere(norm(R(i,:)) - R_m);
    
    ag(i+1,:) = -g * R(i,:)/norm(R(i,:));
    ad(i+1,:) = - (V(i,:) / norm(V(i,:))) * CD * 0.5* rho * norm(V(i,:))^2 * S / m;
    a(i+1,:) = ag(i+1,:) + ad(i+1,:);
    R(i+1,:) = R(i,:) + V(i,:)*dt+a(i,:)*dt^2;
    V(i+1,:) = V(i,:) + a(i,:)*dt;
    
    admag = sqrt(ad(:,1).^2 + ad(:,2).^2 + ad(:,3).^2);
    out.maxaccel = max(admag)/9.81;
        
    % if the s/c has ever been in the atmosphere in_atmus => true
    if (norm(R(i,:)) < (R_m + h_atmos)) 
        out.inatmos = true;
        
    end
        
    % if s/c has been in atmos and is leaving atmos again, stop simulation
    % and check velocity and max deceleration...
    if out.inatmos && ((norm(R(i,:)) > (R_m + h_atmos))) 

        
        %disp(['Vesc = ' num2str(V_esc) ' [m/s], Vend = ' num2str(norm(V(end,:))) ' [m/s]'])
        if (norm(V(end,:)) < V_esc)
            if (out.maxaccel < 3)
               
                out.rightaccel = true;
                
            end
            out.inorbit = true; 
        else
            out.deltav = V_esc-norm(V(end,:));
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


disp(['rx = ' num2str(rx) ' [m], CD = ' num2str(CD) ' [-], in atmosphere: ' num2str(out.inatmos) ', crashed: ' num2str(out.crash) ', in orbit: ' num2str(out.inorbit) ', acceleration: ' num2str(out.maxaccel) ', delta V: ' num2str(out.deltav) ])
end
