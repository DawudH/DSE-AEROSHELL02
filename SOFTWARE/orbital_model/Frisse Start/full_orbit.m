function [ out ] = full_orbit(R0,V0,V_esc,A0,G,M,R_m,h_atm,dt_kep_init,dt_atmos)
%Calculates the full orbit for selected initial conditions until sepcified
%end time

    %Condition to run the script
    orbit = true;

    %Counters
    tp(1) = 0;  %Time
    round = 0; %Number of orbits around mars
    i = 1; %Number of loops
    %%Functions
    [out_hk] = hyperbolic_kepler(R0,V0,A0,G,M,R_m,h_atm,dt_kep_init);
    %%Inputs for while loop
    R(1,:) = out_hk.R;
    speed_sound(1,:) = out_hk.speed_sound;
    V(1,:) = out_hk.V;
    M(1,:) = out_hk.M;
    A(1,:) = out_hk.A;
    Ag(1,:) = out_hk.Ag;
    Ad(1,:) = out_hk.Ad;
    Al(1,:) = out_hk.Al;
    J(1,:) = out_hk.J;
    %Get initial values for conditions
    [out_c] = checks(R,V,V_esc,t);

    %%Functions
    %As long as the s/c is in orbit keep calculating the next position
    while orbit
        %When the s/c is in the atmosphere use the numerical solution including aerodynamic forces
        if out_c.in_atmos

            [out_o] = in_atmosphere( V, R, a, a1, J, atm, CL, CD, dt, R_m, Omega_m, S, m );
            t = t + dt_atmos
        %When the s/c is not in the atmosphere use a kepler orbit
        else
            [out_o] = eliptic_kepler(R,V,A,G,M,R_m,h_atm,dt_kep_init);
            round = round + 1;
            t = t + out_o.t_kep;
        end
        %%New Inputs for while loop
        tp(i+1) = tp(i) + dt_atmos;
        R(i+1,:) = out_o.R;
        speed_sound(i+1,:) = out_o.speed_sound;
        V(i+1,:) = out_o.V;
        M(i+1,:) = out_o.M;
        A(i+1,:) = out_o.A;
        Ag(i+1,:) = out_o.Ag;
        Ad(i+1,:) = out_o.Ad;
        Al(i+1,:) = out_o.Al;
        J(i+1,:) = out_o.J;

        %Function to check when to end the orbit
        [out_c] = checks(R,V,V_esc,t);
        if out_c.crash || out_c.flyby || out_c.t_end
            orbit = false;
        end
        i = i+1;
    end
    
    
    %%Output
    out.R = R;
    out.V = V;
    out.A = A;
    out.Ag = Ag;
    out.Ad = Ad;
    out.Al = Al;
    out.J = J;
end

