function [orbit] = orbit(R,V,a,CD,CL,dt,atm,R_m,Omega_m,S,m)

%Mars atmosphere, get density and g
rho = atm.getDensity(0,0,norm(R)-R_m);
g = atm.getg(norm(R)-R_m);

%gravitational accaleration (vector)
orbit.ag = -g * R/norm(R);
%if in atmosphere there is drag and lift
if rho>0
    %calculate the velocity of the atmosphere
    Vatm = cross(Omega_m,R);
    %calculate the direction of the velocity of the s/c wrt the atmosphere
    vel_unit = (V - Vatm) / norm(V - Vatm);
    %calculate Lift and Drag

    orbit.q = 0.5 * rho * norm(V - Vatm)^2;
    orbit.ad = - vel_unit * CD * orbit.q * S / m;
    orbit.al = cross(vel_unit,[0,0,1]) * CL * orbit.q * S / m;
    
else
    %no lift and drag outside the atmosphere
    orbit.ad = [0,0,0];
    orbit.al = [0,0,0];
    orbit.q = 0;
end
%calculate total accaleration
orbit.a = orbit.ag + orbit.ad + orbit.al;
%calculate new velocity and location
orbit.V = V + a*dt;
orbit.R = R + V*dt + 0.5*a*dt^2;

if (norm(R)-R_m < 100000)
    orbit.speed_sound = atm.getSpeedofsound(0,0, norm(R)-R_m);
    orbit.M = norm(orbit.V) / orbit.speed_sound;
else
    orbit.speed_sound = 0;
    orbit.M = 0;
end
end