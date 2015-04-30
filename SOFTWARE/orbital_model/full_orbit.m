function [t, conditions] = full_orbit(rx,CD)
constants
input_var
conditions.atmos = false; conditions.crash = false; conditions.orbit = true; 
R(1,:) = [rx,ry,0]; %[m]
V(1,:) = [0,-v,0]; %[m/s]
a(1,:) = [0,0,0]; %[m/s^2]
CL = CLCD*CD; %[-]
new = orbit(R,V,a,CD,CL,dt);
V_esc = sqrt(G*M_mars * 2 / (h_atm + R_m)); % escape velocity at the border of the atmosphere
V_esc0 = sqrt(G*M_mars * 2 / norm(R));
i = 1;
h = waitbar(0,'init...');
while true
    %calculate new orbital parameters
    new = orbit(new.R,new.V,new.a,CD,CL,dt);
    R(i+1,:) = new.R;
    V(i+1,:) = new.V;
    a(i+1,:) = new.a;
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
    if (norm(new.R)>(R_m+h_atm)) && (conditions.atmos == true) && (norm(V)>V_esc)
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
    
    if i>t_end/dt
        break;
    end
    i = i+1;
    per = (i*dt/t_end*100);
    waitbar(per/100,h,'%5.1f %%',per);
end
close(h);
t = 0:dt:i*dt;
subplot(3,1,1)
hold on
grid on
plot(t,sqrt(R(:,1).^2 + R(:,2).^2 + R(:,3).^2),'color',cc(1,:));
subplot(3,1,2)
hold on
grid on
plot(t,sqrt(V(:,1).^2 + V(:,2).^2 + V(:,3).^2),'color',cc(1,:));
subplot(3,1,3)
hold on
grid on
plot(t,sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2),'color',cc(1,:));

%circle plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + 104000);
figure('name','Orbit')
grid on
axis equal
hold on
plot(R(:,1),R(:,2))
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g')