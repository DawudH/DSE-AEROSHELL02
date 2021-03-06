function [orbit] = orbit(input,CD,CL,dt,atm,R_m,Omega_m,S,m)
% The input should contain the previous position, velocity, acceleration and jerk
% input.R
% input.V
% input.a
% input.a1 (ai-2)
% input.J

V = input.V;
R = input.R;
a = input.a;
a1 = input.a1;
J = input.J;

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
orbit.J = (a1 - 4*a + 3*orbit.a);

%calculate new velocity and location
orbit.V = V + a*dt + 1/2*J*dt^2;
orbit.R = R + V*dt + 1/2*a*dt^2 + 1/6*J*dt^3;

if (norm(R)-R_m < 100000)
    orbit.speed_sound = atm.getSpeedofsound(0,0, norm(R)-R_m);
    orbit.M = norm(orbit.V) / orbit.speed_sound;
else
    orbit.speed_sound = 0;
    orbit.M = 0;
end
end