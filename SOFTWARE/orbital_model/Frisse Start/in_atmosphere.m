function [ out_o ] = in_atmosphere( V, R, a, a1, J, atm, CL, CD, dt, R_m, omega_m, S, m, phi, Crho )
%IN_ATMOSPHERE Summary of this function goes here
%   Detailed explanation goes here
% a1 = ai-1

%calculate new velocity and location
orbit.V = V + a*dt + 1/2*J*dt^2;
orbit.R = R + V*dt + 1/2*a*dt^2 + 1/6*J*dt^3;

R = orbit.R;
V = orbit.V;

%Mars atmosphere, get density and g
h = norm(R)-R_m;
rho = atm.getCheapDensity(h)*Crho;
g = atm.getg(h);

%gravitational accaleration (vector)
orbit.ag = -g * R/norm(R);
%if in atmosphere there is drag and lift
if rho>0
    %calculate the velocity of the atmosphere
    Vatm = [-omega_m*R(2),omega_m*R(1),0];
    %calculate the direction of the velocity of the s/c wrt the atmosphere
    v_aero = norm(V - Vatm);
    vel_unit = (V - Vatm) / v_aero;
    %calculate Lift and Drag

    orbit.q = 0.5 * rho * v_aero^2;
    orbit.ad = - vel_unit * CD * orbit.q * S / m;
    orbit.aly = CL * orbit.q * S / m*cos(phi);
    orbit.alz = CL * orbit.q * S / m*sin(phi);
    orbit.al = [vel_unit(2),-vel_unit(1),0] * orbit.aly;
    orbit.a_aero = [vel_unit(2),-vel_unit(1),0]* CL * orbit.q * S / m + orbit.ad;
    
    orbit.speed_sound = atm.getCheapSpeedofsound(h);
    orbit.M = v_aero / orbit.speed_sound;
    orbit.T = atm.getCheapTemperature(h);
    
else
    %no lift and drag outside the atmosphere
    orbit.ad = [0,0,0];
    orbit.al = [0,0,0];
    orbit.a_aero = [0,0,0];
    orbit.q = 0;
    
    orbit.speed_sound = 0;
    orbit.M = 0;
    orbit.T = 0;
end


%calculate total accaleration
orbit.a = orbit.ag + orbit.ad + orbit.al;
orbit.J = (a1 - 4*a + 3*orbit.a) / (2*dt);

% output
out_o.speed_sound = orbit.speed_sound;
out_o.R = orbit.R;
out_o.V = orbit.V;
out_o.V_aero = v_aero;
out_o.M = orbit.M;
out_o.A = orbit.a;
out_o.Ag = orbit.ag;
out_o.Ad = orbit.ad;
out_o.Al = orbit.al;
out_o.A_aero = orbit.a_aero;
out_o.J = orbit.J;
out_o.q = orbit.q;
out_o.T = orbit.T;
out_o.rho = rho;
end

