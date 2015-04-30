function [orbit] = orbit(R,V,a,CD,CL,dt,R_m,Omega_m,S,m)

%Mars atmosphere
[ g, p, T, rho, asound ] = mars_atmosphere(norm(R(1,:)) - R_m);
%gravitational accaleration (vector)
orbit.ag = -g * R/norm(R);
%if in atmosphere there is drag and lift
if rho>0
    %calculate the velocity of the atmosphere
    Vatm = cross(Omega_m,R);
    %calculate the direction of the velocity of the s/c wrt the atmosphere
    vel_unit = (V - Vatm) / norm(V - Vatm);
    %calculate Lift and Drag

    orbit.ad = - vel_unit * CD * 0.5* rho * norm(V - Vatm)^2 * S / m;
    orbit.al = - cross(vel_unit,[0,0,1]) * CL * 0.5 * rho * norm(V - Vatm)^2 * S / m;
else
    %no lift and drag outside the atmosphere
    orbit.ad = 0;
    orbit.al = 0;
end
%calculate total accaleration
orbit.a = orbit.ag + orbit.ad + orbit.al;
%calculate new velocity and location
orbit.V = V + a*dt;
orbit.R = R + V*dt + a*dt^2;
end