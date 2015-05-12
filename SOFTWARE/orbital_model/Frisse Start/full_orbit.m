function [ out ] = full_orbit(R0,V0,V_esc,A0,G,M,R_m,h_atm,dt_kep_init,dt_atmos)
%Calculates the full orbit for selected initial conditions until sepcified
%end time
    %%Input
    
    %Condition to run the script
    orbit = true;

    %Counters
    t = 0;  %Time
    round = 0; %Number of orbits around mars

    %%Functions
    [out_hk] = hyperbolic_kepler(R0,V0,A0,G,M,R_m,h_atm,dt_kep_init);
    
    %%Inputs for while loop
    R = out_hk.R;
    V = out_hk.V;
    Vrel = out_hk.Vrel;
    A = out_hk.A;

    %%Functions
    %As long as the s/c is in orbit keep calculating the next position
    while orbit
        %When the s/c is in the atmosphere use the numerical solution including aerodynamic forces
        if out_c.in_atmos
            [out_o] = in_atmosphere( V, R, a, a1, J, atm, CL, CD, dt, R_m, Omega_m, S, m );
            t = t + dt_atmos;
        %When the s/c is not in the atmosphere use a kepler orbit
        else
            [out_o] = eliptic_kepler(R,V,A,G,M,R_m,h_atm,dt_kep_init);
            round = round + 1;
        end
        %%New Inputs for while loop
        R = out_o.R;
        V = out_o.V;
        A = out_o.A;
        J = out_o.J;
        %Function to check when to end the orbit
        [out_c] = checks(R,V,V_esc,t);
        if out_c.crash || out_c.flyby || out_c.t_end
            orbit = false;
        end
    end
    
    
    %%Output
    out.R = 0;
    out.V = 0;
    out.V_rel = 0;
    out.a = 0;
    out.a_g = 0;
    out.a_aero = 0;
    out.hk = out_hk;
end

