function [ out ] = full_orbit(R0,V0,V_esc,A0,G,M,R_m,h_atm,dt_kep_init,dt_atmos)
%Calculates the full orbit for selected initial conditions until sepcified
%end time
    only_hyper = false;
    %%Input
    if only_hyper
        %Condition to run the script
        orbit = true;
        %Note out_c.in_atmos is always true right after the hyperbolic entry is calculated
        out_c.crash = false; out_c.flyby = false; out_c.t_end = false; out_c.in_atmos = true;

        %Counters
        t = 0;  %Time
        round = 0; %Number of orbits around mars
    end
    %%Functions
    [out_hk] = hyperbolic_kepler(R0,V0,A0,G,M,R_m,h_atm,dt_kep_init);
    if only_hyper
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
                [out_o] = in_atmosphere(R,V,Vrel,A,dt_atmos);
                t = t + dt_atmos;
            %When the s/c is not in the atmosphere use a kepler orbit
            else
                [out_o] = eliptic_kepler(R,V,A,G,M,R_m,h_atm,dt_kep_init);
                round = round + 1;
            end
            %%New Inputs for while loop
            R = out_o.R;
            V = out_o.V;
            Vrel = out_o.Vrel;
            A = out_o.A;
            %Function to check when to end the orbit
            [out_c] = checks(R,V,V_esc,t);
            if out_c.crash || out_c.flyby || out_c.t_end
                orbit = false;
            end
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

